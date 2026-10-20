/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useState, useMemo } from 'react';
import { MoleculeRecord, ClusterInfo } from '../types/cheminformatics';
import { Trophy, Download, Search, Check, Sparkles, ArrowUpDown } from 'lucide-react';

interface TopRankedTableProps {
  molecules: MoleculeRecord[];
  clusters: ClusterInfo[];
  activeMoleculeId: string;
  onSelectMolecule: (m: MoleculeRecord) => void;
}

export const TopRankedTable: React.FC<TopRankedTableProps> = ({
  molecules,
  clusters,
  activeMoleculeId,
  onSelectMolecule,
}) => {
  const [searchTerm, setSearchTerm] = useState('');
  const [filterApproved, setFilterApproved] = useState(false);
  const [limit, setLimit] = useState(50);

  // Ranked molecules by SBOF score (descending)
  const sortedMolecules = useMemo(() => {
    return [...molecules].sort((a, b) => {
      const scoreA = a.sbofScore || a.drugScore;
      const scoreB = b.sbofScore || b.drugScore;
      return scoreB - scoreA;
    });
  }, [molecules]);

  const filteredMolecules = useMemo(() => {
    return sortedMolecules
      .filter(m => {
        if (filterApproved && !m.approved && !m.isUserQuery) return false;
        if (!searchTerm) return true;
        const q = searchTerm.toLowerCase();
        return (
          m.name.toLowerCase().includes(q) ||
          m.id.toLowerCase().includes(q) ||
          m.smiles.toLowerCase().includes(q)
        );
      })
      .slice(0, limit);
  }, [sortedMolecules, searchTerm, filterApproved, limit]);

  const exportCSV = () => {
    const headers = ['Rank', 'ID', 'Name', 'Approved', 'SBOF_Score', 'QED', 'MW', 'logP', 'HBA', 'HBD', 'PSA', 'RTB', 'Cluster', 'SMILES'];
    const rows = sortedMolecules.map((m, idx) => [
      idx + 1,
      m.id,
      `"${m.name.replace(/"/g, '""')}"`,
      m.approved ? 'TRUE' : 'FALSE',
      (m.sbofScore || m.drugScore).toFixed(4),
      m.descriptors.qed,
      m.descriptors.mw,
      m.descriptors.logP,
      m.descriptors.hba,
      m.descriptors.hbd,
      m.descriptors.psa,
      m.descriptors.rtb,
      m.cluster + 1,
      `"${m.smiles}"`,
    ]);
    const csvContent = [headers.join(','), ...rows.map(r => r.join(','))].join('\n');
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', `sbof_top_ranked_molecules_${Date.now()}.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  return (
    <div className="bg-slate-900 border border-slate-800 rounded-lg flex flex-col font-mono text-xs overflow-hidden h-full">
      {/* Table Header & Controls */}
      <div className="p-3 border-b border-slate-800 bg-slate-950/60 flex flex-col sm:flex-row items-center justify-between gap-2.5">
        <div className="flex items-center gap-2">
          <Trophy className="w-4 h-4 text-amber-400" />
          <span className="font-semibold text-slate-200">SBOF Re-Ranked Chemical Library</span>
          <span className="text-[11px] text-slate-500">({molecules.length} entities)</span>
        </div>

        <div className="flex items-center gap-2 w-full sm:w-auto">
          {/* Search box */}
          <div className="relative flex-1 sm:w-44">
            <Search className="w-3.5 h-3.5 text-slate-500 absolute left-2.5 top-2" />
            <input
              type="text"
              placeholder="Search drug / ID..."
              value={searchTerm}
              onChange={e => setSearchTerm(e.target.value)}
              className="w-full bg-slate-900 border border-slate-700/80 rounded pl-7 pr-2 py-1 text-slate-200 text-[11px] focus:outline-none focus:border-cyan-500"
            />
          </div>

          {/* FDA Filter */}
          <button
            onClick={() => setFilterApproved(!filterApproved)}
            className={`px-2 py-1 text-[11px] rounded border transition ${
              filterApproved
                ? 'bg-emerald-500/20 text-emerald-300 border-emerald-500/50'
                : 'bg-slate-800 text-slate-400 border-slate-700 hover:text-slate-200'
            }`}
          >
            FDA Only
          </button>

          {/* Export CSV */}
          <button
            onClick={exportCSV}
            className="flex items-center gap-1 px-2.5 py-1 bg-cyan-950/80 hover:bg-cyan-900 text-cyan-300 border border-cyan-800/80 rounded transition text-[11px]"
            title="Download CSV for computational biology pipeline"
          >
            <Download className="w-3 h-3" />
            <span>CSV</span>
          </button>
        </div>
      </div>

      {/* Table Content */}
      <div className="flex-1 overflow-auto">
        <table className="w-full text-left border-collapse">
          <thead className="bg-slate-950/80 sticky top-0 z-10 border-b border-slate-800 text-[11px] text-slate-400">
            <tr>
              <th className="py-2 px-3 w-12 text-center">Rank</th>
              <th className="py-2 px-3">Compound / ID</th>
              <th className="py-2 px-3 text-right">SBOF Score</th>
              <th className="py-2 px-3 text-right">QED</th>
              <th className="py-2 px-3 text-right">MW (Da)</th>
              <th className="py-2 px-3 text-right">logP</th>
              <th className="py-2 px-3 text-center">Cluster</th>
              <th className="py-2 px-3 text-center">Status</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-800/60 text-slate-300">
            {filteredMolecules.map((mol, idx) => {
              const isActive = mol.id === activeMoleculeId;
              const isUser = mol.isUserQuery;
              const cluster = clusters.find(c => c.id === mol.cluster);

              return (
                <tr
                  key={mol.id}
                  onClick={() => onSelectMolecule(mol)}
                  className={`cursor-pointer transition-colors ${
                    isActive
                      ? 'bg-cyan-950/60 border-l-4 border-l-cyan-400 text-white'
                      : isUser
                      ? 'bg-indigo-950/40 hover:bg-indigo-900/40 border-l-2 border-l-indigo-400'
                      : 'hover:bg-slate-800/50'
                  }`}
                >
                  <td className="py-1.5 px-3 text-center font-bold">
                    {idx === 0 && <span className="text-amber-400">🥇 1</span>}
                    {idx === 1 && <span className="text-slate-300">🥈 2</span>}
                    {idx === 2 && <span className="text-amber-600">🥉 3</span>}
                    {idx > 2 && <span className="text-slate-500">#{idx + 1}</span>}
                  </td>
                  <td className="py-1.5 px-3">
                    <div className="font-semibold flex items-center gap-1.5 text-slate-200">
                      {isUser && <Sparkles className="w-3 h-3 text-cyan-400 animate-pulse" />}
                      <span>{mol.name}</span>
                    </div>
                    <div className="text-[10px] text-slate-500">{mol.id}</div>
                  </td>
                  <td className="py-1.5 px-3 text-right font-bold text-cyan-300">
                    {(mol.sbofScore || mol.drugScore).toFixed(4)}
                  </td>
                  <td className="py-1.5 px-3 text-right text-slate-300">
                    {mol.descriptors.qed}
                  </td>
                  <td className="py-1.5 px-3 text-right text-slate-300">
                    {mol.descriptors.mw}
                  </td>
                  <td className="py-1.5 px-3 text-right text-slate-300">
                    {mol.descriptors.logP}
                  </td>
                  <td className="py-1.5 px-3 text-center">
                    <span
                      className="px-1.5 py-0.5 rounded text-[10px] font-medium"
                      style={{
                        backgroundColor: `${cluster?.color || '#06b6d4'}25`,
                        color: cluster?.color || '#06b6d4',
                        border: `1px solid ${cluster?.color || '#06b6d4'}50`,
                      }}
                    >
                      C{mol.cluster + 1}
                    </span>
                  </td>
                  <td className="py-1.5 px-3 text-center">
                    {mol.approved ? (
                      <span className="text-[10px] px-1.5 py-0.2 bg-emerald-950/80 text-emerald-400 border border-emerald-800/80 rounded">
                        Approved
                      </span>
                    ) : (
                      <span className="text-[10px] text-slate-500">Investigational</span>
                    )}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
};
