import { useState, useRef, useEffect, useCallback } from "react";

/**
 * FHN Trajectory Animator — modular interactive panel
 * ---------------------------------------------------
 * Full 2D stochastic FitzHugh–Nagumo, integrated with Euler–Maruyama,
 * matching kernel.py exactly:
 *
 *   dv = (v - v^3/3 - w + I) dt + sigma * sqrt(dt) * xi
 *   dw = eps * (v + a - b w) dt
 *
 * Live phase plane (with nullclines + fixed points + stability readout)
 * and a scrolling voltage trace. Sliders for I, eps, sigma, dt, speed.
 *
 * Defaults: a = 0.7, b = 0.8 (project standard). Drop-in: no props required.
 */

// ----------------------------------------------------------------------------
// Model helpers (pure JS mirror of kernel.py geometry)
// ----------------------------------------------------------------------------

const cubicNull = (v, I) => v - (v * v * v) / 3 - 0 + I; // w = v - v^3/3 + I  (v-nullcline)
const lineNull = (v, a, b) => (v + a) / b; //                                  (w-nullcline)

// g(v) = 0 at fixed points: v - v^3/3 - (v+a)/b + I
const gFP = (v, I, a, b) => v - (v * v * v) / 3 - (v + a) / b + I;

function findFixedPoints(I, a, b) {
  const pts = [];
  const lo = -3,
    hi = 3,
    N = 2000;
  let vPrev = lo,
    gPrev = gFP(lo, I, a, b);
  for (let k = 1; k <= N; k++) {
    const v = lo + ((hi - lo) * k) / N;
    const g = gFP(v, I, a, b);
    if (gPrev === 0) pts.push(vPrev);
    else if (gPrev * g < 0) {
      // bisection
      let a0 = vPrev,
        b0 = v,
        fa = gPrev;
      for (let it = 0; it < 60; it++) {
        const m = 0.5 * (a0 + b0);
        const fm = gFP(m, I, a, b);
        if (fa * fm <= 0) b0 = m;
        else {
          a0 = m;
          fa = fm;
        }
      }
      pts.push(0.5 * (a0 + b0));
    }
    vPrev = v;
    gPrev = g;
  }
  return pts;
}

function classifyFP(v, eps, b) {
  // J = [[1 - v^2, -1], [eps, -eps*b]]
  const tr = 1 - v * v - eps * b;
  const det = eps * (1 - b + b * v * v); // (1-v^2)(-eps b) + eps
  const disc = tr * tr - 4 * det;
  let kind;
  if (det < 0) kind = "saddle";
  else {
    const osc = disc < 0;
    if (tr < 0) kind = osc ? "stable spiral" : "stable node";
    else if (tr > 0) kind = osc ? "unstable spiral" : "unstable node";
    else kind = "center (Hopf)";
  }
  return { tr, det, disc, kind };
}

// Regime label from the stable/unstable structure of the unique relevant FP.
function regimeOf(I, eps, a, b) {
  const fps = findFixedPoints(I, a, b);
  // pick the lowest-v FP as the "operating" one (continuation from excitable)
  const classified = fps.map((v) => ({ v, ...classifyFP(v, eps, b) }));
  const anyUnstable = classified.some(
    (c) => c.det > 0 && c.tr > 0 && c.v > -1 && c.v < 1
  );
  // single-FP cases
  if (anyUnstable) return { label: "Tonic (limit cycle)", fps: classified };
  // operating FP
  const op = classified[0];
  if (!op) return { label: "—", fps: classified };
  if (op.v < -1) return { label: "Excitable (left branch)", fps: classified };
  if (op.v > 1) return { label: "Excitable (right branch)", fps: classified };
  if (op.kind.startsWith("stable spiral"))
    return { label: "Resonator (stable spiral)", fps: classified };
  return { label: op.kind, fps: classified };
}

// ----------------------------------------------------------------------------
// Gaussian noise (Box–Muller)
// ----------------------------------------------------------------------------
let _spare = null;
function gauss() {
  if (_spare !== null) {
    const s = _spare;
    _spare = null;
    return s;
  }
  let u = 0,
    v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  const mag = Math.sqrt(-2 * Math.log(u));
  _spare = mag * Math.sin(2 * Math.PI * v);
  return mag * Math.cos(2 * Math.PI * v);
}

// ----------------------------------------------------------------------------
// Parameter control: type-in number box + continuous slider.
// Module-level (NOT redefined each render) so the text input keeps focus.
// The slider is bounded for convenient dragging; the number box accepts ANY
// value (including outside the slider range) so you can explore freely.
// ----------------------------------------------------------------------------
function fmtNum(v) {
  if (!isFinite(v)) return String(v);
  // tidy display without clobbering the stored full-precision value
  return String(Number(v.toPrecision(5)));
}

function ParamControl({ label, value, min, max, step, onChange }) {
  const [text, setText] = useState(fmtNum(value));
  useEffect(() => {
    if (parseFloat(text) !== value) setText(fmtNum(value));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [value]);

  const commit = (raw) => {
    const v = parseFloat(raw);
    if (!isNaN(v)) onChange(v);
  };

  return (
    <div className="flex flex-col gap-1">
      <div className="flex justify-between items-center text-xs text-slate-300">
        <span>{label}</span>
        <input
          type="number"
          value={text}
          step={step}
          onChange={(e) => {
            setText(e.target.value);
            commit(e.target.value);
          }}
          className="w-24 px-1.5 py-0.5 text-right font-mono text-slate-100 bg-slate-800 border border-slate-700 rounded focus:outline-none focus:border-emerald-400"
        />
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step="any"
        value={Math.min(max, Math.max(min, value))}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="w-full accent-emerald-400"
      />
    </div>
  );
}

// ----------------------------------------------------------------------------
// Component
// ----------------------------------------------------------------------------
export default function FHNAnimator() {
  const [I, setI] = useState(-0.1);
  const [eps, setEps] = useState(0.08);
  const [sigma, setSigma] = useState(0.1);
  const [dt, setDt] = useState(0.005);
  const [stepsPerFrame, setSteps] = useState(40);
  const [running, setRunning] = useState(true);
  const [showNull, setShowNull] = useState(true);

  const [a, setA] = useState(0.7);
  const [b, setB] = useState(0.8);
  // custom trajectory start point (so reset doesn't auto-snap to a stable FP)
  const [v0, setV0] = useState(-1.2);
  const [w0, setW0] = useState(-0.6);

  // mutable simulation state (avoid re-render per step)
  const state = useRef({ v: -1.2, w: -0.6, t: 0, trail: [], trace: [], spikes: 0, lastAbove: false });
  const phaseRef = useRef(null);
  const traceRef = useRef(null);
  const rafRef = useRef(null);
  const paramsRef = useRef({ I, eps, sigma, dt, stepsPerFrame, a, b });
  paramsRef.current = { I, eps, sigma, dt, stepsPerFrame, a, b };

  const reg = regimeOf(I, eps, a, b);
  // ---- live bifurcation diagnostics ----
  const fpList = reg.fps; // [{v, tr, det, kind}]
  const vH = 1 - eps * b > 0 ? Math.sqrt(1 - eps * b) : null; // Hopf threshold |v|
  const vSN = b > 1 ? Math.sqrt(1 - 1 / b) : null; //            saddle-node |v| (needs b>1)
  const epsB2 = eps * b * b; //                                  εb²: >1 clean fold, <1 BT tangle
  const dotColor = (kind) =>
    kind.startsWith("stable") ? "#36d399" : kind === "saddle" ? "#fbbd23" : "#f87272";

  // reset from the user-specified (v0, w0) instead of a stable fixed point,
  // so the trajectory isn't immediately absorbed.
  const resetSim = useCallback(() => {
    state.current = {
      v: v0,
      w: w0,
      t: 0,
      trail: [],
      trace: [],
      spikes: 0,
      lastAbove: false,
    };
  }, [v0, w0]);

  // snap the start point to the lowest-v fixed point on its w-nullcline
  const startFromFP = useCallback(() => {
    const fps = findFixedPoints(I, a, b);
    const vv = fps.length ? fps[0] : -1.2;
    setV0(vv);
    setW0((vv + a) / b);
  }, [I, a, b]);

  // ---- phase-plane drawing helpers ----
  const VMIN = -2.6,
    VMAX = 2.6,
    WMIN = -1.4,
    WMAX = 2.4;

  // click anywhere on the phase plane to set the start point there
  const onPhaseClick = useCallback(
    (e) => {
      const cv = phaseRef.current;
      if (!cv) return;
      const rect = cv.getBoundingClientRect();
      const px = ((e.clientX - rect.left) / rect.width) * cv.width;
      const py = ((e.clientY - rect.top) / rect.height) * cv.height;
      const v = VMIN + (px / cv.width) * (VMAX - VMIN);
      const w = WMIN + (1 - py / cv.height) * (WMAX - WMIN);
      setV0(v);
      setW0(w);
      // place trajectory there immediately
      state.current = { v, w, t: 0, trail: [], trace: [], spikes: 0, lastAbove: false };
    },
    []
  );

  function drawPhase() {
    const cv = phaseRef.current;
    if (!cv) return;
    const ctx = cv.getContext("2d");
    const W = cv.width,
      H = cv.height;
    const x = (v) => ((v - VMIN) / (VMAX - VMIN)) * W;
    const y = (w) => H - ((w - WMIN) / (WMAX - WMIN)) * H;

    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = "#0b1020";
    ctx.fillRect(0, 0, W, H);

    // grid
    ctx.strokeStyle = "rgba(255,255,255,0.06)";
    ctx.lineWidth = 1;
    for (let gv = -2; gv <= 2; gv++) {
      ctx.beginPath();
      ctx.moveTo(x(gv), 0);
      ctx.lineTo(x(gv), H);
      ctx.stroke();
    }
    for (let gw = -1; gw <= 2; gw++) {
      ctx.beginPath();
      ctx.moveTo(0, y(gw));
      ctx.lineTo(W, y(gw));
      ctx.stroke();
    }

    const { I, a, b, eps } = paramsRef.current;
    if (showNull) {
      // v-nullcline (cubic)
      ctx.strokeStyle = "#4ea1ff";
      ctx.lineWidth = 2;
      ctx.beginPath();
      for (let i = 0; i <= 300; i++) {
        const v = VMIN + ((VMAX - VMIN) * i) / 300;
        const w = cubicNull(v, I);
        if (i === 0) ctx.moveTo(x(v), y(w));
        else ctx.lineTo(x(v), y(w));
      }
      ctx.stroke();
      // w-nullcline (line)
      ctx.strokeStyle = "#ff7a59";
      ctx.beginPath();
      for (let i = 0; i <= 2; i++) {
        const v = i === 0 ? VMIN : VMAX;
        const w = lineNull(v, a, b);
        if (i === 0) ctx.moveTo(x(v), y(w));
        else ctx.lineTo(x(v), y(w));
      }
      ctx.stroke();
    }

    // trail
    const tr = state.current.trail;
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = "rgba(120,220,160,0.85)";
    ctx.beginPath();
    for (let i = 0; i < tr.length; i++) {
      const p = tr[i];
      if (i === 0) ctx.moveTo(x(p[0]), y(p[1]));
      else ctx.lineTo(x(p[0]), y(p[1]));
    }
    ctx.stroke();

    // fixed points
    findFixedPoints(I, a, b).forEach((v) => {
      const c = classifyFP(v, paramsRef.current.eps, b);
      const w = lineNull(v, a, b);
      ctx.beginPath();
      ctx.arc(x(v), y(w), 5, 0, 2 * Math.PI);
      ctx.fillStyle = c.kind.startsWith("stable")
        ? "#36d399"
        : c.kind === "saddle"
        ? "#fbbd23"
        : "#f87272";
      ctx.fill();
      ctx.strokeStyle = "#0b1020";
      ctx.stroke();
    });

    // chosen start point
    ctx.beginPath();
    ctx.arc(x(v0), y(w0), 6.5, 0, 2 * Math.PI);
    ctx.strokeStyle = "rgba(255,255,255,0.85)";
    ctx.lineWidth = 1.5;
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x(v0) - 4, y(w0));
    ctx.lineTo(x(v0) + 4, y(w0));
    ctx.moveTo(x(v0), y(w0) - 4);
    ctx.lineTo(x(v0), y(w0) + 4);
    ctx.stroke();

    // current point
    const { v, w } = state.current;
    ctx.beginPath();
    ctx.arc(x(v), y(w), 4.5, 0, 2 * Math.PI);
    ctx.fillStyle = "#ffffff";
    ctx.fill();

    // axes labels
    ctx.fillStyle = "rgba(255,255,255,0.55)";
    ctx.font = "12px ui-sans-serif, system-ui";
    ctx.fillText("v", W - 16, H - 8);
    ctx.fillText("w", 8, 14);
  }

  function drawTrace() {
    const cv = traceRef.current;
    if (!cv) return;
    const ctx = cv.getContext("2d");
    const W = cv.width,
      H = cv.height;
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = "#0b1020";
    ctx.fillRect(0, 0, W, H);

    const data = state.current.trace;
    const TMIN = -2.6,
      TMAX = 2.6;
    const y = (val) => H - ((val - TMIN) / (TMAX - TMIN)) * H;

    // threshold line (v = 1)
    ctx.strokeStyle = "rgba(255,255,255,0.18)";
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(0, y(1));
    ctx.lineTo(W, y(1));
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.strokeStyle = "#7ee0a0";
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    const n = data.length;
    for (let i = 0; i < n; i++) {
      const px = (i / Math.max(1, n - 1)) * W;
      if (i === 0) ctx.moveTo(px, y(data[i]));
      else ctx.lineTo(px, y(data[i]));
    }
    ctx.stroke();

    ctx.fillStyle = "rgba(255,255,255,0.55)";
    ctx.font = "12px ui-sans-serif, system-ui";
    ctx.fillText("v(t)", 8, 14);
  }

  // ---- animation loop ----
  useEffect(() => {
    function frame() {
      const { I, eps, sigma, dt, stepsPerFrame, a, b } = paramsRef.current;
      const s = state.current;
      const sq = Math.sqrt(dt);
      for (let k = 0; k < stepsPerFrame; k++) {
        const v = s.v,
          w = s.w;
        const dv = (v - (v * v * v) / 3 - w + I) * dt + sigma * sq * gauss();
        const dw = eps * (v + a - b * w) * dt;
        s.v = v + dv;
        s.w = w + dw;
        s.t += dt;
        // spike detection (upward cross of v = 1)
        const above = s.v >= 1;
        if (above && !s.lastAbove) s.spikes += 1;
        s.lastAbove = above;
        // record
        s.trail.push([s.v, s.w]);
        if (s.trail.length > 1200) s.trail.shift();
        s.trace.push(s.v);
        if (s.trace.length > 1000) s.trace.shift();
      }
      drawPhase();
      drawTrace();
      rafRef.current = requestAnimationFrame(frame);
    }
    if (running) rafRef.current = requestAnimationFrame(frame);
    return () => cancelAnimationFrame(rafRef.current);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [running, showNull, a, b]);

  // redraw static parts when params change while paused
  useEffect(() => {
    if (!running) {
      drawPhase();
      drawTrace();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [I, eps, sigma, dt, showNull, running, a, b, v0, w0]);

  return (
    <div className="w-full max-w-5xl mx-auto p-4 bg-slate-950 rounded-2xl text-slate-100">
      <div className="flex items-center justify-between mb-3">
        <h2 className="text-lg font-semibold">FHN Trajectory Animator</h2>
        <span className="text-xs px-2 py-1 rounded-full bg-slate-800 border border-slate-700">
          regime:&nbsp;
          <span className="font-medium text-emerald-300">{reg.label}</span>
        </span>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        <div className="rounded-xl overflow-hidden border border-slate-800 relative">
          <canvas
            ref={phaseRef}
            width={460}
            height={380}
            onClick={onPhaseClick}
            title="Click to set the starting point (v0, w0)"
            className="w-full block cursor-crosshair"
          />
        </div>
        <div className="flex flex-col gap-4">
          <div className="rounded-xl overflow-hidden border border-slate-800">
            <canvas ref={traceRef} width={460} height={170} className="w-full block" />
          </div>

          <div className="grid grid-cols-2 gap-x-4 gap-y-3">
            <ParamControl label="I (current)" value={I} min={-0.5} max={1.6} step={0.005} onChange={setI} />
            <ParamControl label="ε (timescale)" value={eps} min={0.005} max={0.5} step={0.005} onChange={setEps} />
            <ParamControl label="σ (noise)" value={sigma} min={0} max={0.5} step={0.005} onChange={setSigma} />
            <ParamControl label="dt" value={dt} min={0.001} max={0.02} step={0.001} onChange={setDt} />
            <ParamControl label="a (recovery)" value={a} min={0} max={1.5} step={0.01} onChange={setA} />
            <ParamControl label="b (recovery)" value={b} min={0.05} max={3.5} step={0.01} onChange={setB} />
            <ParamControl label="v₀ (start v)" value={v0} min={VMIN} max={VMAX} step={0.01} onChange={setV0} />
            <ParamControl label="w₀ (start w)" value={w0} min={WMIN} max={WMAX} step={0.01} onChange={setW0} />
            <ParamControl label="speed (steps/frame)" value={stepsPerFrame} min={1} max={120} step={1} onChange={(x) => setSteps(Math.round(x))} />
          </div>

          <div className="flex items-center gap-2 flex-wrap">
            <button
              onClick={() => setRunning((r) => !r)}
              className="px-3 py-1.5 rounded-lg bg-emerald-500 hover:bg-emerald-400 text-slate-950 text-sm font-medium"
            >
              {running ? "Pause" : "Play"}
            </button>
            <button
              onClick={resetSim}
              className="px-3 py-1.5 rounded-lg bg-slate-700 hover:bg-slate-600 text-sm"
            >
              Reset
            </button>
            <button
              onClick={() => {
                startFromFP();
                const fps = findFixedPoints(I, a, b);
                const vv = fps.length ? fps[0] : -1.2;
                state.current = {
                  v: vv,
                  w: (vv + a) / b,
                  t: 0,
                  trail: [],
                  trace: [],
                  spikes: 0,
                  lastAbove: false,
                };
              }}
              className="px-3 py-1.5 rounded-lg bg-slate-700 hover:bg-slate-600 text-sm"
            >
              Start at FP
            </button>
            <button
              onClick={() => setShowNull((s) => !s)}
              className="px-3 py-1.5 rounded-lg bg-slate-700 hover:bg-slate-600 text-sm"
            >
              {showNull ? "Hide" : "Show"} nullclines
            </button>
            <span className="text-xs text-slate-400 ml-auto font-mono">
              spikes: {state.current.spikes}
            </span>
          </div>
        </div>
      </div>

      <div className="mt-3 grid grid-cols-1 md:grid-cols-2 gap-3 text-xs">
        <div className="rounded-xl border border-slate-800 p-3 bg-slate-900/40">
          <div className="text-slate-300 font-medium mb-1">
            Fixed points ({fpList.length})
          </div>
          <div className="flex flex-col gap-0.5 font-mono">
            {fpList.map((f, i) => (
              <div key={i} className="flex items-center gap-2">
                <span
                  className="inline-block w-2.5 h-2.5 rounded-full shrink-0"
                  style={{ background: dotColor(f.kind) }}
                />
                <span className="text-slate-200 w-20">v = {f.v.toFixed(3)}</span>
                <span className="text-slate-400">{f.kind}</span>
              </div>
            ))}
          </div>
        </div>
        <div className="rounded-xl border border-slate-800 p-3 bg-slate-900/40 font-mono text-slate-300 flex flex-col gap-1">
          <div>
            regime: <span className="text-emerald-300">{reg.label}</span>
          </div>
          <div>
            Hopf threshold |v_H| = {vH !== null ? vH.toFixed(3) : "—"}
            <span className="text-slate-500"> (stable if |v*| {">"} |v_H|)</span>
          </div>
          <div>
            saddle-node |v_sn| ={" "}
            {vSN !== null ? vSN.toFixed(3) : `— (need b ${">"} 1)`}
          </div>
          <div className="mt-1">
            εb² = {epsB2.toFixed(3)} →{" "}
            <span className={epsB2 > 1 ? "text-emerald-300" : "text-amber-300"}>
              {epsB2 > 1
                ? "clean bistable fold (well stable up to SN)"
                : "Bogdanov–Takens tangle (Hopf pre-empts SN)"}
            </span>
          </div>
        </div>
      </div>

      <div className="mt-3 text-[11px] text-slate-400 leading-relaxed">
        Blue = v-nullcline (cubic), orange = w-nullcline (line). White crosshair = chosen start point; click the phase plane or edit v₀/w₀ to set it. Dots: green = stable,
        red = unstable, yellow = saddle. Euler–Maruyama on
        dv = (v − v³/3 − w + I)dt + σ√dt·ξ, dw = ε(v + a − bw)dt, with a = {a}, b = {b}.
      </div>
    </div>
  );
}
