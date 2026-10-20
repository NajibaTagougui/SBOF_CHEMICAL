/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React from 'react';
import { MolecularDescriptors, DrugScoreBreakdown } from '../types/cheminformatics';
import { DESCRIPTOR_BOUNDS } from '../services/sbofEngine';
import { CheckCircle2, AlertTriangle, XCircle, ShieldCheck, Zap, Activity } from 'lucide-react';

interface DescriptorPanelProps {
  descriptors: MolecularDescriptors;
  drugScore: DrugScoreBreakdown;
  compoundName: string;
}

export const DescriptorPanel: React.FC<DescriptorPanelProps> = ({
  descriptors,
  drugScore,
  compoundName,
}) => {
  const getDescriptorPercent = (val: number, min: number, max: number) => {
    return Math.max(0, Math.min(100, ((val - min) / (max - min)) * 100));
  };

  const descriptorItems = [
    {
      key: 'MW',
      label: 'Molecular Weight',
      value: `${descriptors.mw} Da`,
      raw: descriptors.mw,
      range: '100 – 894 Da',
      percent: getDescriptorPercent(descriptors.mw, 100, 900),
      optimal: descriptors.mw <= 500,
      ruleText: 'Lipinski ≤ 500',
    },
    {
      key: 'logP',
      label: 'Partition Coeff (logP)',
      value: `${descriptors.logP}`,
      raw: descriptors.logP,
      range: '-2.0 – 6.0',
      percent: getDescriptorPercent(descriptors.logP, -4, 8),
      optimal: descriptors.logP <= 5.0 && descriptors.logP >= -1.0,
      ruleText: 'Lipinski ≤ 5.0',
    },
    {
      key: 'QED',
      label: 'Quantitative Drug-Likeness',
      value: `${descriptors.qed}`,
      raw: descriptors.qed,
      range: '0.0 – 1.0',
      percent: descriptors.qed * 100,
      optimal: descriptors.qed >= 0.5,
      ruleText: 'Bickerton Score',
    },
    {
      key: 'HBA',
      label: 'H-Bond Acceptors',
      value: `${descriptors.hba}`,
      raw: descriptors.hba,
      range: '0 – 15',
      percent: getDescriptorPercent(descriptors.hba, 0, 16),
      optimal: descriptors.hba <= 10,
      ruleText: 'Lipinski ≤ 10',
    },
    {
      key: 'HBD',
      label: 'H-Bond Donors',
      value: `${descriptors.hbd}`,
      raw: descriptors.hbd,
      range: '0 – 8',
      percent: getDescriptorPercent(descriptors.hbd, 0, 10),
      optimal: descriptors.hbd <= 5,
      ruleText: 'Lipinski ≤ 5',
    },
    {
      key: 'PSA',
      label: 'Polar Surface Area (TPSA)',
      value: `${descriptors.psa} Å²`,
      raw: descriptors.psa,
      range: '0 – 200 Å²',
      percent: getDescriptorPercent(descriptors.psa, 0, 220),
      optimal: descriptors.psa <= 140,
      ruleText: 'Veber ≤ 140 Å²',
    },
    {
      key: 'RTB',
      label: 'Rotatable Bonds',
      value: `${descriptors.rtb}`,
      raw: descriptors.rtb,
      range: '0 – 15',
      percent: getDescriptorPercent(descriptors.rtb, 0, 20),
      optimal: descriptors.rtb <= 10,
      ruleText: 'Veber ≤ 10',
    },
  ];

  return (
    <div className="bg-slate-900 border border-slate-800 rounded-lg p-3.5 flex flex-col gap-3 font-mono text-xs">
      {/* Header: Score Summary */}
      <div className="flex items-center justify-between border-b border-slate-800 pb-2.5">
        <div className="flex items-center gap-2">
          <Activity className="w-4 h-4 text-cyan-400" />
          <span className="font-semibold text-slate-200">7 RDKit-Equivalent Descriptors</span>
        </div>
        <div className="flex items-center gap-2">
          <span className="text-[11px] text-slate-400">SBOF Rank:</span>
          <span className="px-2 py-0.5 bg-cyan-950 text-cyan-300 border border-cyan-800 rounded font-bold">
            #{drugScore.sbofRank} (Top {drugScore.percentile}%)
          </span>
        </div>
      </div>

      {/* SBOF Score Gauge */}
      <div className="p-2.5 bg-slate-950/70 border border-slate-800/80 rounded flex items-center justify-between">
        <div>
          <div className="text-[10px] text-slate-400 uppercase tracking-wide">SBOF Learned Drug-Likeness Potential</div>
          <div className="text-xl font-bold text-transparent bg-clip-text bg-gradient-to-r from-cyan-400 to-emerald-400">
            {drugScore.sbofScore.toFixed(4)}
          </div>
        </div>
        <div className="flex gap-2 text-right">
          <div>
            <div className="text-[10px] text-slate-400">Lipinski Ro5</div>
            <div className={`flex items-center gap-1 font-semibold ${drugScore.lipinskiPassed ? 'text-emerald-400' : 'text-amber-400'}`}>
              {drugScore.lipinskiPassed ? <CheckCircle2 className="w-3.5 h-3.5" /> : <AlertTriangle className="w-3.5 h-3.5" />}
              <span>{drugScore.lipinskiPassed ? 'Passed' : `${drugScore.lipinskiViolations} Viol.`}</span>
            </div>
          </div>
          <div className="border-l border-slate-800 pl-2">
            <div className="text-[10px] text-slate-400">Veber Bioavail.</div>
            <div className={`flex items-center gap-1 font-semibold ${drugScore.veberPassed ? 'text-emerald-400' : 'text-amber-400'}`}>
              {drugScore.veberPassed ? <CheckCircle2 className="w-3.5 h-3.5" /> : <AlertTriangle className="w-3.5 h-3.5" />}
              <span>{drugScore.veberPassed ? 'Passed' : `${drugScore.veberViolations} Viol.`}</span>
            </div>
          </div>
        </div>
      </div>

      {/* Descriptors Metric Bars */}
      <div className="space-y-2">
        {descriptorItems.map(item => (
          <div key={item.key} className="space-y-1">
            <div className="flex items-center justify-between text-[11px]">
              <span className="text-slate-300 font-medium">
                {item.key} <span className="text-slate-500 text-[10px]">({item.label})</span>
              </span>
              <div className="flex items-center gap-2">
                <span className={`text-[10px] px-1.5 py-0.2 rounded border ${
                  item.optimal
                    ? 'bg-emerald-950/60 text-emerald-400 border-emerald-900/60'
                    : 'bg-amber-950/60 text-amber-400 border-amber-900/60'
                }`}>
                  {item.ruleText}
                </span>
                <span className="font-bold text-white min-w-[55px] text-right">{item.value}</span>
              </div>
            </div>
            {/* Progress Bar */}
            <div className="w-full h-1.5 bg-slate-800 rounded-full overflow-hidden">
              <div
                className={`h-full transition-all duration-300 rounded-full ${
                  item.optimal ? 'bg-gradient-to-r from-cyan-500 to-emerald-400' : 'bg-gradient-to-r from-amber-500 to-red-400'
                }`}
                style={{ width: `${Math.max(4, Math.min(100, item.percent))}%` }}
              />
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};
