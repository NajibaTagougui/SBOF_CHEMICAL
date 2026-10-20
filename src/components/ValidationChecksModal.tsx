/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useState } from 'react';
import { runDatasetValidationChecks, SBOF_BENCHMARKS } from '../services/sbofEngine';
import { MoleculeRecord } from '../types/cheminformatics';
import { ShieldCheck, CheckCircle2, AlertTriangle, XCircle, TrendingUp, X } from 'lucide-react';

interface ValidationChecksModalProps {
  isOpen: boolean;
  onClose: () => void;
  dataset: MoleculeRecord[];
}

export const ValidationChecksModal: React.FC<ValidationChecksModalProps> = ({
  isOpen,
  onClose,
  dataset,
}) => {
  const [activeTab, setActiveTab] = useState<'validator' | 'benchmarks'>('validator');

  if (!isOpen) return null;

  const checks = runDatasetValidationChecks(dataset);
  const passCount = checks.filter(c => c.status === 'Pass').length;
  const warnCount = checks.filter(c => c.status === 'Warning').length;
  const expectedFailCount = checks.filter(c => c.status.includes('Fail')).length;

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/70 backdrop-blur-sm p-4 animate-in fade-in duration-200 font-mono text-xs">
      <div className="bg-slate-900 border border-slate-700/80 rounded-xl shadow-2xl max-w-4xl w-full max-h-[85vh] flex flex-col overflow-hidden">
        {/* Modal Header */}
        <div className="p-4 border-b border-slate-800 bg-slate-950 flex items-center justify-between">
          <div className="flex items-center gap-2.5">
            <ShieldCheck className="w-5 h-5 text-emerald-400" />
            <div>
              <h2 className="text-sm font-bold text-slate-100">
                DrugBank 5.1.10 Quality Validation & Benchmark Suite
              </h2>
              <p className="text-[11px] text-slate-400">
                27-Point Dataset Validator (Table 3) & SBOF Baseline Comparisons (Table 4)
              </p>
            </div>
          </div>
          <button
            onClick={onClose}
            className="p-1.5 text-slate-400 hover:text-white rounded hover:bg-slate-800 transition"
          >
            <X className="w-4 h-4" />
          </button>
        </div>

        {/* Tab Controls */}
        <div className="flex items-center justify-between px-4 py-2 border-b border-slate-800 bg-slate-900/90">
          <div className="flex gap-2">
            <button
              onClick={() => setActiveTab('validator')}
              className={`px-3 py-1.5 rounded transition flex items-center gap-1.5 text-[11px] ${
                activeTab === 'validator'
                  ? 'bg-emerald-500/20 text-emerald-300 border border-emerald-500/50 font-bold'
                  : 'text-slate-400 hover:text-white'
              }`}
            >
              <CheckCircle2 className="w-3.5 h-3.5" />
              <span>27-Point Validator ({passCount}/27 Passed)</span>
            </button>
            <button
              onClick={() => setActiveTab('benchmarks')}
              className={`px-3 py-1.5 rounded transition flex items-center gap-1.5 text-[11px] ${
                activeTab === 'benchmarks'
                  ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/50 font-bold'
                  : 'text-slate-400 hover:text-white'
              }`}
            >
              <TrendingUp className="w-3.5 h-3.5" />
              <span>Baseline Comparisons (Table 4)</span>
            </button>
          </div>

          <div className="flex items-center gap-3 text-[11px] text-slate-400">
            <span>Pass: <strong className="text-emerald-400">{passCount}</strong></span>
            <span>Warn: <strong className="text-amber-400">{warnCount}</strong></span>
            <span>Expected Fail: <strong className="text-slate-300">{expectedFailCount}</strong></span>
          </div>
        </div>

        {/* Body */}
        <div className="flex-1 p-4 overflow-auto bg-slate-950">
          {activeTab === 'validator' ? (
            <div className="space-y-1.5">
              <div className="text-[11px] text-slate-400 mb-2 italic">
                Checks reproduce the exact 27 criteria from Tagougui & Ben Abdallah (Table 3), verifying absence of NaNs, molecular weight bounded support, descriptor variability, Lipinski/Veber rule compliance, and rank-7 non-degeneracy.
              </div>
              <div className="border border-slate-800 rounded overflow-hidden">
                <table className="w-full text-left border-collapse">
                  <thead className="bg-slate-900 border-b border-slate-800 text-[11px] text-slate-400">
                    <tr>
                      <th className="py-2 px-3">Validation Metric</th>
                      <th className="py-2 px-3">Value / Result</th>
                      <th className="py-2 px-3 text-center">Status</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-800/60 text-slate-300">
                    {checks.map((c, i) => (
                      <tr key={i} className="hover:bg-slate-900/40">
                        <td className="py-1.5 px-3 font-medium text-slate-200">{c.name}</td>
                        <td className="py-1.5 px-3 text-slate-400">{c.value}</td>
                        <td className="py-1.5 px-3 text-center">
                          {c.status === 'Pass' && (
                            <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-[10px] bg-emerald-950 text-emerald-400 border border-emerald-800">
                              <CheckCircle2 className="w-3 h-3" /> Pass
                            </span>
                          )}
                          {c.status === 'Warning' && (
                            <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-[10px] bg-amber-950 text-amber-400 border border-amber-800">
                              <AlertTriangle className="w-3 h-3" /> Warning
                            </span>
                          )}
                          {c.status.includes('Fail') && (
                            <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-[10px] bg-slate-800 text-slate-300 border border-slate-700">
                              Expected ({c.status.split(' ')[0]})
                            </span>
                          )}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          ) : (
            <div className="space-y-4">
              <div className="p-3 bg-cyan-950/40 border border-cyan-800/60 rounded text-cyan-200 text-[11px] leading-relaxed">
                <strong>Algorithmic Superiority of SBOF:</strong> Even against an optimally tuned Sinkhorn baseline (ϵ = 0.005, TV = 0.1753), SBOF retains a <strong>91.9% structural improvement</strong> in Total Variation distance (TV = 0.0141, Pearson r = 0.9996), avoiding the n × n memory bottleneck through single-potential mirror descent.
              </div>

              <div className="border border-slate-800 rounded overflow-hidden">
                <table className="w-full text-left border-collapse">
                  <thead className="bg-slate-900 border-b border-slate-800 text-[11px] text-slate-400">
                    <tr>
                      <th className="py-2 px-3">Method</th>
                      <th className="py-2 px-3 text-right">TV Distance ↓</th>
                      <th className="py-2 px-3 text-right">KL Divergence ↓</th>
                      <th className="py-2 px-3 text-right">Pearson r ↑</th>
                      <th className="py-2 px-3 text-right">Runtime (s)</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-800/60 text-slate-300">
                    {SBOF_BENCHMARKS.map((b, i) => (
                      <tr
                        key={i}
                        className={b.highlight ? 'bg-cyan-950/40 font-semibold text-cyan-300 border-l-4 border-l-cyan-400' : 'hover:bg-slate-900/40'}
                      >
                        <td className="py-2 px-3 flex items-center gap-1.5">
                          {b.highlight && <ShieldCheck className="w-3.5 h-3.5 text-cyan-400" />}
                          <span>{b.method}</span>
                        </td>
                        <td className="py-2 px-3 text-right">{b.tv.toFixed(4)}</td>
                        <td className="py-2 px-3 text-right">{b.kl}</td>
                        <td className="py-2 px-3 text-right">{b.pearson.toFixed(4)}</td>
                        <td className="py-2 px-3 text-right text-slate-400">{b.runtime}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-3 text-[11px]">
                <div className="p-3 bg-slate-900/60 border border-slate-800 rounded">
                  <div className="font-bold text-slate-200 mb-1">Mirror Descent Update (Eq 4)</div>
                  <code className="text-cyan-300 block bg-slate-950 p-2 rounded mb-2">
                    g_{`{k+1}`} = g_k - η · (log m_k - log ν)_centered, η = 0.03
                  </code>
                  <p className="text-slate-400">
                    Avoids the n × n transport matrix coupling P, directly updating the scalar dual potential per molecule on the probability simplex.
                  </p>
                </div>

                <div className="p-3 bg-slate-900/60 border border-slate-800 rounded">
                  <div className="font-bold text-slate-200 mb-1">Target Tempered Softmax (Eq 2)</div>
                  <code className="text-emerald-300 block bg-slate-950 p-2 rounded mb-2">
                    ν_i = exp(β · s_i) / Σ exp(β · s_j), β = 4.0
                  </code>
                  <p className="text-slate-400">
                    Shapes the target drug-likeness measure from composite score s_i, interpolating between uniform and sharp drug-like modes.
                  </p>
                </div>
              </div>
            </div>
          )}
        </div>

        {/* Footer */}
        <div className="p-3 border-t border-slate-800 bg-slate-950 text-slate-400 text-[11px] flex items-center justify-between">
          <div>Reference: Tagougui & Ben Abdallah, 2026 (DrugBank 5.1.10)</div>
          <button
            onClick={onClose}
            className="px-3 py-1 bg-slate-800 hover:bg-slate-700 text-slate-200 rounded transition"
          >
            Close
          </button>
        </div>
      </div>
    </div>
  );
};
