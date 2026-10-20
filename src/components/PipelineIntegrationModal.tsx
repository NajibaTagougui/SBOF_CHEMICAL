/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useState } from 'react';
import { PipelineExportPayload, Conformer3D } from '../types/cheminformatics';
import { Terminal, Copy, Check, Download, Code, FileText, X } from 'lucide-react';

interface PipelineIntegrationModalProps {
  isOpen: boolean;
  onClose: () => void;
  payload: PipelineExportPayload;
  conformer: Conformer3D | null;
}

export const PipelineIntegrationModal: React.FC<PipelineIntegrationModalProps> = ({
  isOpen,
  onClose,
  payload,
  conformer,
}) => {
  const [activeTab, setActiveTab] = useState<'json' | 'python' | 'snakemake' | 'curl' | 'sdf'>('json');
  const [copied, setCopied] = useState(false);

  if (!isOpen) return null;

  const jsonString = JSON.stringify(payload, null, 2);

  // Generate SDF representation from conformer coordinates
  const generateSDFString = () => {
    if (!conformer || conformer.atoms.length === 0) return 'No 3D conformer loaded.';
    const header = `${payload.queryCompound.name}\n  SBOF-Engine-v1.0\n3D Conformer Coordinates\n`;
    const numAtoms = String(conformer.atoms.length).padStart(3, ' ');
    const numBonds = String(conformer.bonds.length).padStart(3, ' ');
    const countsLine = `${numAtoms}${numBonds}  0  0  0  0  0  0  0  0999 V2000\n`;

    let atomLines = '';
    for (const a of conformer.atoms) {
      const x = a.x.toFixed(4).padStart(10, ' ');
      const y = a.y.toFixed(4).padStart(10, ' ');
      const z = a.z.toFixed(4).padStart(10, ' ');
      const sym = a.symbol.padEnd(3, ' ');
      atomLines += `${x}${y}${z} ${sym} 0  0  0  0  0  0  0  0  0\n`;
    }

    let bondLines = '';
    for (const b of conformer.bonds) {
      const src = String(b.source + 1).padStart(3, ' ');
      const tgt = String(b.target + 1).padStart(3, ' ');
      const ord = String(b.order === 1.5 ? 4 : Math.round(b.order)).padStart(3, ' ');
      bondLines += `${src}${tgt}${ord}  0  0  0  0\n`;
    }

    const properties = `M  END\n> <SBOF_POTENTIAL>\n${payload.queryCompound.sbofMetrics.sbofPotential}\n\n> <QED>\n${payload.queryCompound.descriptors.qed}\n\n> <PCA_COORDINATES>\n${payload.queryCompound.pcaCoordinates.pc1},${payload.queryCompound.pcaCoordinates.pc2},${payload.queryCompound.pcaCoordinates.pc3}\n\n$$$$\n`;

    return header + countsLine + atomLines + bondLines + properties;
  };

  const sdfString = generateSDFString();

  // Python integration script
  const pythonSnippet = `#!/usr/bin/env python3
"""
Computational Biology Pipeline Integration: SBOF Drug-Likeness Engine
Project: Schrödinger Bridge Onsager Flow on DrugBank 5.1.10
"""

import json
from dataclasses import dataclass
from typing import Dict, Any

@dataclass
class SBOFPayload:
    smiles: str
    descriptors: Dict[str, float]
    pca_coordinates: Dict[str, float]
    cluster_id: int
    sbof_potential: float
    rank: int

def parse_sbof_output(raw_json: str) -> SBOFPayload:
    data = json.loads(raw_json)["queryCompound"]
    return SBOFPayload(
        smiles=data["smiles"],
        descriptors=data["descriptors"],
        pca_coordinates=data["pcaCoordinates"],
        cluster_id=data["clusterAssignment"]["clusterId"],
        sbof_potential=data["sbofMetrics"]["sbofPotential"],
        rank=data["sbofMetrics"]["libraryRank"]
    )

# Direct execution hook in pipeline:
if __name__ == "__main__":
    query_smiles = "${payload.queryCompound.smiles}"
    print(f"[Pipeline] Processing compound: {query_smiles}")
    print(f"[Pipeline] SBOF Potential: {${payload.queryCompound.sbofMetrics.sbofPotential}} (Rank #{${payload.queryCompound.sbofMetrics.libraryRank}})")
    print(f"[Pipeline] Cluster Assignment: Cluster {${payload.queryCompound.clusterAssignment.clusterId + 1}}")
    print(f"[Pipeline] Drug-Likeness Lipinski Passed: {${payload.queryCompound.druglikenessRules.lipinskiRuleOfFive.passed}}")
`;

  // Snakemake workflow rule
  const snakemakeSnippet = `rule sbof_chemical_space_clustering:
    input:
        smiles = "data/raw_compounds/{sample}.smi",
        drugbank_space = "models/sbof_drugbank_5.1.10.h5"
    output:
        descriptors = "results/{sample}/descriptors.json",
        conformer = "results/{sample}/conformer_3d.sdf",
        pca_coords = "results/{sample}/pca_projection.tsv",
        sbof_ranking = "results/{sample}/top_ranked.tsv"
    params:
        sigma = 0.35,
        beta = 4.0,
        eta = 0.03
    conda:
        "envs/rdkit_sbof.yaml"
    shell:
        """
        sbof-cli cluster \\
            --input {input.smiles} \\
            --db-reference {input.drugbank_space} \\
            --sigma {params.sigma} \\
            --beta {params.beta} \\
            --output-json {output.descriptors} \\
            --output-sdf {output.conformer}
        """
`;

  // cURL endpoint
  const curlSnippet = `# Direct CLI / API invocation for automated pipelines
curl -X POST "https://api.sbof-drugbank.org/v1/cluster-rank" \\
  -H "Content-Type: application/json" \\
  -d '{
    "smiles": "${payload.queryCompound.smiles}",
    "descriptors": {
      "mw": ${payload.queryCompound.descriptors.mw},
      "logP": ${payload.queryCompound.descriptors.logP},
      "qed": ${payload.queryCompound.descriptors.qed},
      "hba": ${payload.queryCompound.descriptors.hba},
      "hbd": ${payload.queryCompound.descriptors.hbd},
      "psa": ${payload.queryCompound.descriptors.psa},
      "rtb": ${payload.queryCompound.descriptors.rtb}
    },
    "k_clusters": 2,
    "beta": 4.0,
    "sigma": 0.35
  }'
`;

  const getCurrentText = () => {
    switch (activeTab) {
      case 'json': return jsonString;
      case 'python': return pythonSnippet;
      case 'snakemake': return snakemakeSnippet;
      case 'curl': return curlSnippet;
      case 'sdf': return sdfString;
    }
  };

  const copyToClipboard = () => {
    navigator.clipboard.writeText(getCurrentText());
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const downloadFile = () => {
    const isSdf = activeTab === 'sdf';
    const isJson = activeTab === 'json';
    const filename = isSdf
      ? `${payload.queryCompound.name.toLowerCase().replace(/\s+/g, '_')}_conformer.sdf`
      : isJson
      ? `${payload.queryCompound.name.toLowerCase().replace(/\s+/g, '_')}_sbof_payload.json`
      : `sbof_${activeTab}_script.txt`;

    const blob = new Blob([getCurrentText()], { type: isJson ? 'application/json' : 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/70 backdrop-blur-sm p-4 animate-in fade-in duration-200 font-mono text-xs">
      <div className="bg-slate-900 border border-slate-700/80 rounded-xl shadow-2xl max-w-4xl w-full max-h-[85vh] flex flex-col overflow-hidden">
        {/* Modal Header */}
        <div className="p-4 border-b border-slate-800 bg-slate-950 flex items-center justify-between">
          <div className="flex items-center gap-2.5">
            <Terminal className="w-5 h-5 text-cyan-400" />
            <div>
              <h2 className="text-sm font-bold text-slate-100">Direct Computational Pipeline Integration</h2>
              <p className="text-[11px] text-slate-400">
                Direct export, API payload, and workflow scripts for virtual screening pipelines
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

        {/* Tab Bar */}
        <div className="flex items-center justify-between px-4 py-2 border-b border-slate-800 bg-slate-900/90 flex-wrap gap-2">
          <div className="flex gap-1.5">
            <button
              onClick={() => setActiveTab('json')}
              className={`px-3 py-1.5 rounded transition flex items-center gap-1.5 text-[11px] ${
                activeTab === 'json' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/50 font-bold' : 'text-slate-400 hover:text-white'
              }`}
            >
              <FileText className="w-3.5 h-3.5" />
              <span>JSON Payload</span>
            </button>
            <button
              onClick={() => setActiveTab('sdf')}
              className={`px-3 py-1.5 rounded transition flex items-center gap-1.5 text-[11px] ${
                activeTab === 'sdf' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/50 font-bold' : 'text-slate-400 hover:text-white'
              }`}
            >
              <Code className="w-3.5 h-3.5" />
              <span>3D SDF Conformer</span>
            </button>
            <button
              onClick={() => setActiveTab('python')}
              className={`px-3 py-1.5 rounded transition flex items-center gap-1.5 text-[11px] ${
                activeTab === 'python' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/50 font-bold' : 'text-slate-400 hover:text-white'
              }`}
            >
              <Terminal className="w-3.5 h-3.5" />
              <span>Python Script</span>
            </button>
            <button
              onClick={() => setActiveTab('snakemake')}
              className={`px-3 py-1.5 rounded transition flex items-center gap-1.5 text-[11px] ${
                activeTab === 'snakemake' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/50 font-bold' : 'text-slate-400 hover:text-white'
              }`}
            >
              <Code className="w-3.5 h-3.5" />
              <span>Snakemake / Nextflow</span>
            </button>
            <button
              onClick={() => setActiveTab('curl')}
              className={`px-3 py-1.5 rounded transition flex items-center gap-1.5 text-[11px] ${
                activeTab === 'curl' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/50 font-bold' : 'text-slate-400 hover:text-white'
              }`}
            >
              <Terminal className="w-3.5 h-3.5" />
              <span>cURL REST</span>
            </button>
          </div>

          <div className="flex items-center gap-2">
            <button
              onClick={copyToClipboard}
              className="flex items-center gap-1.5 px-3 py-1.5 bg-slate-800 hover:bg-slate-700 text-slate-200 rounded border border-slate-700 text-[11px] transition"
            >
              {copied ? <Check className="w-3.5 h-3.5 text-emerald-400" /> : <Copy className="w-3.5 h-3.5" />}
              <span>{copied ? 'Copied' : 'Copy'}</span>
            </button>
            <button
              onClick={downloadFile}
              className="flex items-center gap-1.5 px-3 py-1.5 bg-cyan-600 hover:bg-cyan-500 text-slate-950 font-bold rounded shadow transition text-[11px]"
            >
              <Download className="w-3.5 h-3.5" />
              <span>Download</span>
            </button>
          </div>
        </div>

        {/* Code Content Viewport */}
        <div className="flex-1 p-4 overflow-auto bg-slate-950">
          <pre className="text-slate-300 font-mono text-[11px] leading-relaxed whitespace-pre select-all">
            {getCurrentText()}
          </pre>
        </div>

        {/* Footer */}
        <div className="p-3 border-t border-slate-800 bg-slate-950 text-slate-400 text-[11px] flex items-center justify-between">
          <div>
            Target: <span className="text-slate-200 font-semibold">{payload.queryCompound.name}</span> |
            SBOF Potential: <span className="text-cyan-400 font-semibold">{payload.queryCompound.sbofMetrics.sbofPotential}</span>
          </div>
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
