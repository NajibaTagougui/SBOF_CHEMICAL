/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useState, useEffect, useMemo, useRef } from 'react';
import {
  MoleculeRecord,
  MolecularDescriptors,
  Conformer3D,
  ClusterInfo,
  DrugScoreBreakdown,
  PipelineExportPayload,
} from './types/cheminformatics';
import {
  calculateDescriptors,
  evaluateDrugLikeness,
  generateConformerFromSmiles,
  parseCoordinateFile,
} from './services/cheminformatics';
import {
  DRUGBANK_LIBRARY,
  normalizeDescriptors,
  projectPCA,
  evaluateSBOFNet,
  solveSBOF,
  clusterChemicalSpace,
  SBOF_CONFIG,
} from './services/sbofEngine';
import { PCAViewer3D } from './components/PCAViewer3D';
import { Conformer3DViewer } from './components/Conformer3DViewer';
import { DescriptorPanel } from './components/DescriptorPanel';
import { TopRankedTable } from './components/TopRankedTable';
import { PipelineIntegrationModal } from './components/PipelineIntegrationModal';
import { ValidationChecksModal } from './components/ValidationChecksModal';
import {
  Atom,
  Binary,
  Cpu,
  Download,
  FileCode,
  FileUp,
  FlaskConical,
  Layers,
  Play,
  RotateCw,
  Search,
  ShieldCheck,
  Sparkles,
  Terminal,
  Upload,
} from 'lucide-react';

// Preset molecules from paper
const PRESET_MOLECULES = [
  { name: 'Levomenthol', smiles: 'CC1CCC(C(C1)O)C(C)C', note: 'Top SBOF Paper Hit' },
  { name: 'Menadione (Vit K3)', smiles: 'CC1=C(C(=O)c2ccccc2C1=O)C', note: 'High SBOF Potential' },
  { name: 'Thymol', smiles: 'Cc1ccc(c(c1)O)C(C)C', note: 'Natural Phenol Lead' },
  { name: 'Ciprofloxacin', smiles: 'C1CC1N2C=C(C(=O)c3cc(c(cc32)N4CCNCC4)F)C(=O)O', note: 'Antibiotic Drug' },
  { name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O', note: 'Classic NSAID' },
  { name: 'Imatinib', smiles: 'Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C', note: 'Kinase Inhibitor' },
  { name: 'Artemisinin', smiles: 'CC1CCC2C(C(=O)OC3(C2(C1)C(=O)O3)C)C', note: 'Antimalarial Lead' },
  { name: 'Remdesivir', smiles: 'CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)c2ccc3n2ncnc3N)O)O)Oc4ccccc4', note: 'High-MW Complex' },
];

export default function App() {
  // Input states
  const [smilesInput, setSmilesInput] = useState('CC1CCC(C(C1)O)C(C)C'); // Levomenthol default
  const [compoundName, setCompoundName] = useState('Levomenthol');
  const [fileInputActive, setFileInputActive] = useState(false);
  const [pastedCoords, setPastedCoords] = useState('');
  const [activeTab3D, setActiveTab3D] = useState<'both' | 'pca' | 'conformer'>('both');

  // Clustering & SBOF states
  const [kClusters, setKClusters] = useState(2); // optimal k=2 from paper
  const [betaSoftmax, setBetaSoftmax] = useState(4.0); // tempered softmax beta=4.0
  const [library, setLibrary] = useState<MoleculeRecord[]>(DRUGBANK_LIBRARY);
  const [clusterData, setClusterData] = useState<{
    clusters: ClusterInfo[];
    silhouetteScore: number;
    daviesBouldinIndex: number;
  }>({ clusters: [], silhouetteScore: 0.479, daviesBouldinIndex: 0.933 });

  // Active compound states
  const [activeDescriptors, setActiveDescriptors] = useState<MolecularDescriptors>({
    mw: 156.3, logP: 2.85, qed: 0.684, hba: 1, hbd: 1, psa: 20.2, rtb: 1
  });
  const [activeConformer, setActiveConformer] = useState<Conformer3D | null>(null);
  const [activeRecord, setActiveRecord] = useState<MoleculeRecord | null>(null);
  const [isSolving, setIsSolving] = useState(false);
  const [solvingMessage, setSolvingMessage] = useState('');

  // Modals
  const [isPipelineModalOpen, setIsPipelineModalOpen] = useState(false);
  const [isValidationsModalOpen, setIsValidationsModalOpen] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Initialize Library & Initial Clustering on mount
  useEffect(() => {
    runClusteringAndScoring(DRUGBANK_LIBRARY, kClusters, betaSoftmax);
    processMolecularInput('CC1CCC(C(C1)O)C(C)C', 'Levomenthol');
  }, []);

  // Update clusters when k changes
  const handleKChange = (newK: number) => {
    setKClusters(newK);
    const updated = [...library];
    const clustered = clusterChemicalSpace(updated, newK);
    setClusterData(clustered);
    setLibrary(updated);
  };

  // Full SBOF Dual-Ascent & Clustering calculation
  const runClusteringAndScoring = (baseLibrary: MoleculeRecord[], k: number, beta: number) => {
    setIsSolving(true);
    setSolvingMessage('Computing Gaussian Similarity Kernel K (σ=0.35)...');

    setTimeout(() => {
      setSolvingMessage('Running Centered Mirror-Descent Dual Ascent (η=0.03)...');
      const sbofResult = solveSBOF(baseLibrary, beta, 50);

      // Assign scores to library
      const updatedLib = baseLibrary.map((m, idx) => ({
        ...m,
        sbofScore: sbofResult.scores[idx] || m.drugScore,
      }));

      setSolvingMessage('Projecting 3D PCA chemical space & K-Means clusters...');
      const clustered = clusterChemicalSpace(updatedLib, k);

      setLibrary(updatedLib);
      setClusterData(clustered);
      setIsSolving(false);
      setSolvingMessage('');
    }, 150);
  };

  // Process user input (SMILES or coordinate-derived)
  const processMolecularInput = (smiles: string, name: string, preParsedConformer?: Conformer3D) => {
    const descriptors = calculateDescriptors(smiles, preParsedConformer);
    const conformer = preParsedConformer || generateConformerFromSmiles(smiles);
    const norm = normalizeDescriptors(descriptors);
    const pca = projectPCA(norm);
    const drugEval = evaluateDrugLikeness(descriptors);
    const sbofNetScore = evaluateSBOFNet(descriptors);

    setActiveDescriptors(descriptors);
    setActiveConformer(conformer);
    setCompoundName(name);

    // Assign to nearest cluster
    let bestCluster = 0;
    let bestDist = Infinity;
    if (clusterData.clusters.length > 0) {
      clusterData.clusters.forEach(c => {
        const dx = pca[0] - c.centroid[0];
        const dy = pca[1] - c.centroid[1];
        const dz = pca[2] - c.centroid[2];
        const d = dx * dx + dy * dy + dz * dz;
        if (d < bestDist) {
          bestDist = d;
          bestCluster = c.id;
        }
      });
    }

    const userRecord: MoleculeRecord = {
      id: `QUERY-${Math.floor(1000 + Math.random() * 9000)}`,
      name,
      smiles,
      approved: false,
      descriptors,
      pca,
      cluster: bestCluster,
      sbofScore: sbofNetScore,
      drugScore: drugEval.compositeScore,
      isUserQuery: true,
    };

    setActiveRecord(userRecord);

    // Check if compound already in library or merge it
    setLibrary(prev => {
      const filtered = prev.filter(m => !m.isUserQuery);
      return [userRecord, ...filtered];
    });
  };

  // Find top nearest neighbors in chemical space
  const nearestNeighbors = useMemo(() => {
    if (!activeRecord) return [];
    const [ax, ay, az] = activeRecord.pca;
    return library
      .filter(m => m.id !== activeRecord.id)
      .map(m => {
        const dx = m.pca[0] - ax;
        const dy = m.pca[1] - ay;
        const dz = m.pca[2] - az;
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        return {
          id: m.id,
          name: m.name,
          distance: Number(dist.toFixed(3)),
          sbofScore: m.sbofScore || m.drugScore,
          approved: m.approved,
        };
      })
      .sort((a, b) => a.distance - b.distance)
      .slice(0, 5);
  }, [activeRecord, library]);

  // Drug Score Breakdown
  const drugScoreBreakdown: DrugScoreBreakdown = useMemo(() => {
    const evalResult = evaluateDrugLikeness(activeDescriptors);
    const sorted = [...library].sort((a, b) => (b.sbofScore || b.drugScore) - (a.sbofScore || a.drugScore));
    const rank = activeRecord ? sorted.findIndex(m => m.id === activeRecord.id) + 1 : 1;
    const percentile = Number(((rank / Math.max(1, library.length)) * 100).toFixed(1));

    return {
      compositeScore: evalResult.compositeScore,
      lipinskiPassed: evalResult.lipinskiPassed,
      lipinskiViolations: evalResult.lipinskiViolations,
      veberPassed: evalResult.veberPassed,
      veberViolations: evalResult.veberViolations,
      sbofScore: activeRecord?.sbofScore || evaluateSBOFNet(activeDescriptors),
      sbofRank: rank > 0 ? rank : 1,
      percentile,
    };
  }, [activeDescriptors, activeRecord, library]);

  // Pipeline Export Payload
  const pipelinePayload: PipelineExportPayload = useMemo(() => {
    return {
      version: '1.0.0',
      timestamp: new Date().toISOString(),
      engine: 'Schrodinger-Bridge-Onsager-Flow (SBOF)',
      source: 'DrugBank 5.1.10 Reference Chemical Space',
      queryCompound: {
        name: compoundName,
        smiles: smilesInput,
        descriptors: activeDescriptors,
        pcaCoordinates: {
          pc1: activeRecord?.pca[0] || 0,
          pc2: activeRecord?.pca[1] || 0,
          pc3: activeRecord?.pca[2] || 0,
        },
        clusterAssignment: {
          clusterId: activeRecord?.cluster || 0,
          clusterName: clusterData.clusters.find(c => c.id === activeRecord?.cluster)?.name || 'Cluster 1',
          silhouetteOptimalK: kClusters,
        },
        sbofMetrics: {
          compositeScore: drugScoreBreakdown.compositeScore,
          sbofPotential: drugScoreBreakdown.sbofScore,
          libraryRank: drugScoreBreakdown.sbofRank,
          totalLibrarySize: library.length,
          percentile: drugScoreBreakdown.percentile,
        },
        druglikenessRules: {
          lipinskiRuleOfFive: {
            passed: drugScoreBreakdown.lipinskiPassed,
            violations: drugScoreBreakdown.lipinskiViolations,
          },
          veberOralBioavailability: {
            passed: drugScoreBreakdown.veberPassed,
            violations: drugScoreBreakdown.veberViolations,
          },
        },
        nearestNeighbors: nearestNeighbors.map(nn => ({
          id: nn.id,
          name: nn.name,
          distance: nn.distance,
          sbofScore: nn.sbofScore,
          approved: nn.approved,
        })),
      },
    };
  }, [compoundName, smilesInput, activeDescriptors, activeRecord, clusterData, kClusters, drugScoreBreakdown, library.length, nearestNeighbors]);

  // Handle File Upload
  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = event => {
      const content = event.target?.result as string;
      if (content) {
        try {
          const { conformer, derivedSmiles } = parseCoordinateFile(content, file.name);
          const cleanName = file.name.replace(/\.[^/.]+$/, '');
          setSmilesInput(derivedSmiles);
          processMolecularInput(derivedSmiles, cleanName, conformer);
        } catch (err) {
          console.error('File parsing error:', err);
        }
      }
    };
    reader.readAsText(file);
    e.target.value = '';
  };

  // Handle Pasted Coordinates
  const handlePastedCoordinatesSubmit = () => {
    if (!pastedCoords.trim()) return;
    try {
      const { conformer, derivedSmiles } = parseCoordinateFile(pastedCoords, 'pasted.sdf');
      setSmilesInput(derivedSmiles);
      processMolecularInput(derivedSmiles, 'Pasted-Conformer', conformer);
      setFileInputActive(false);
      setPastedCoords('');
    } catch (err) {
      console.error('Failed to parse pasted coordinates', err);
    }
  };

  // Handle selecting any molecule in the 3D space or table
  const handleSelectMolecule = (m: MoleculeRecord) => {
    setSmilesInput(m.smiles);
    setCompoundName(m.name);
    setActiveDescriptors(m.descriptors);
    const conf = generateConformerFromSmiles(m.smiles);
    setActiveConformer(conf);
    setActiveRecord(m);
  };

  return (
    <div className="flex flex-col h-screen w-screen bg-slate-950 text-slate-100 font-mono select-none overflow-hidden">
      {/* Top Header */}
      <header className="h-13 bg-slate-900 border-b border-slate-800 px-4 flex items-center justify-between z-20 shrink-0">
        <div className="flex items-center gap-3">
          <div className="w-8 h-8 rounded bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center shadow-lg shadow-cyan-500/20">
            <Atom className="w-5 h-5 text-slate-950 font-bold" />
          </div>
          <div>
            <div className="flex items-center gap-2">
              <h1 className="text-sm font-bold tracking-tight text-white">
                SBOF Chemical Space Clustering & Top Drug-Likeness Engine
              </h1>
              <span className="px-1.5 py-0.2 text-[10px] bg-cyan-950 text-cyan-300 border border-cyan-800 rounded font-semibold">
                DrugBank 5.1.10
              </span>
            </div>
            <p className="text-[10px] text-slate-400 hidden md:block">
              Schrödinger Bridge Onsager Flow • 7 RDKit Descriptors • 3D Atomic Conformer • 3D PCA Space • k-Means
            </p>
          </div>
        </div>

        {/* Global Action Buttons */}
        <div className="flex items-center gap-2">
          {/* Re-solve SBOF */}
          <button
            onClick={() => runClusteringAndScoring(library, kClusters, betaSoftmax)}
            disabled={isSolving}
            className="flex items-center gap-1.5 px-2.5 py-1.5 bg-slate-800 hover:bg-slate-700 text-slate-200 border border-slate-700 rounded text-xs transition"
            title="Re-run SBOF Centered Mirror Descent"
          >
            <RotateCw className={`w-3.5 h-3.5 text-cyan-400 ${isSolving ? 'animate-spin' : ''}`} />
            <span className="hidden sm:inline">Solve SBOF</span>
          </button>

          {/* 27-Point Validation */}
          <button
            onClick={() => setIsValidationsModalOpen(true)}
            className="flex items-center gap-1.5 px-2.5 py-1.5 bg-emerald-950/70 hover:bg-emerald-900 text-emerald-300 border border-emerald-800/80 rounded text-xs transition"
          >
            <ShieldCheck className="w-3.5 h-3.5 text-emerald-400" />
            <span className="hidden sm:inline">27-Point Checks</span>
          </button>

          {/* Direct Computational Pipeline Integration */}
          <button
            onClick={() => setIsPipelineModalOpen(true)}
            className="flex items-center gap-1.5 px-3 py-1.5 bg-gradient-to-r from-cyan-600 to-blue-600 hover:from-cyan-500 hover:to-blue-500 text-slate-950 font-bold rounded shadow-md text-xs transition"
          >
            <Terminal className="w-3.5 h-3.5" />
            <span>Pipeline Integration</span>
          </button>
        </div>
      </header>

      {/* Input Ribbon */}
      <div className="bg-slate-900/90 border-b border-slate-800/80 px-4 py-2 flex flex-col gap-2 shrink-0 z-20">
        <div className="flex flex-col md:flex-row items-center gap-2">
          {/* Preset Buttons */}
          <div className="flex items-center gap-1.5 overflow-x-auto w-full md:w-auto shrink-0 pb-1 md:pb-0">
            <span className="text-[11px] text-slate-400 font-semibold uppercase tracking-wider flex items-center gap-1">
              <Sparkles className="w-3 h-3 text-cyan-400" /> Examples:
            </span>
            {PRESET_MOLECULES.map(p => (
              <button
                key={p.name}
                onClick={() => {
                  setSmilesInput(p.smiles);
                  processMolecularInput(p.smiles, p.name);
                }}
                className={`px-2 py-0.5 text-[11px] rounded transition shrink-0 border ${
                  compoundName === p.name
                    ? 'bg-cyan-500 text-slate-950 font-bold border-cyan-400'
                    : 'bg-slate-800 text-slate-300 border-slate-700 hover:bg-slate-700'
                }`}
                title={p.note}
              >
                {p.name}
              </button>
            ))}
          </div>

          {/* SMILES Input Bar */}
          <div className="flex-1 flex items-center gap-2 w-full">
            <div className="relative flex-1">
              <input
                type="text"
                value={smilesInput}
                onChange={e => setSmilesInput(e.target.value)}
                onKeyDown={e => {
                  if (e.key === 'Enter') {
                    processMolecularInput(smilesInput, 'Custom Compound');
                  }
                }}
                placeholder="Enter SMILES (e.g. CC1CCC(C(C1)O)C(C)C)..."
                className="w-full bg-slate-950 border border-slate-700/80 rounded px-3 py-1.5 text-xs text-cyan-300 focus:outline-none focus:border-cyan-400 font-mono"
              />
            </div>

            <button
              onClick={() => processMolecularInput(smilesInput, 'Custom Compound')}
              className="px-3 py-1.5 bg-cyan-600 hover:bg-cyan-500 text-slate-950 font-bold text-xs rounded transition flex items-center gap-1.5 shrink-0"
            >
              <Play className="w-3.5 h-3.5 fill-current" />
              <span>Calculate & Cluster</span>
            </button>

            {/* Upload File Button */}
            <input
              type="file"
              ref={fileInputRef}
              onChange={handleFileUpload}
              accept=".sdf,.pdb,.mol,.xyz"
              className="hidden"
            />
            <button
              onClick={() => fileInputRef.current?.click()}
              className="flex items-center gap-1 px-2.5 py-1.5 bg-slate-800 hover:bg-slate-700 text-slate-200 border border-slate-700 rounded text-xs transition shrink-0"
              title="Upload 3D Structural File (SDF, PDB, MOL, XYZ)"
            >
              <FileUp className="w-3.5 h-3.5 text-emerald-400" />
              <span className="hidden sm:inline">Upload 3D</span>
            </button>

            {/* Paste Coords Toggle */}
            <button
              onClick={() => setFileInputActive(!fileInputActive)}
              className={`flex items-center gap-1 px-2.5 py-1.5 text-xs rounded border transition shrink-0 ${
                fileInputActive
                  ? 'bg-indigo-500/20 text-indigo-300 border-indigo-500/50'
                  : 'bg-slate-800 text-slate-400 border-slate-700 hover:text-slate-200'
              }`}
              title="Paste 3D structural coordinates"
            >
              <FileCode className="w-3.5 h-3.5" />
              <span className="hidden sm:inline">Paste Coords</span>
            </button>
          </div>
        </div>

        {/* Paste Box Drawer if active */}
        {fileInputActive && (
          <div className="p-3 bg-slate-950 border border-slate-800 rounded flex flex-col gap-2 animate-in fade-in">
            <div className="flex items-center justify-between text-[11px] text-slate-400">
              <span>Paste 3D coordinates (SDF, PDB, MOL, or XYZ format):</span>
              <button
                onClick={() => setFileInputActive(false)}
                className="text-slate-500 hover:text-slate-300"
              >
                Cancel
              </button>
            </div>
            <textarea
              value={pastedCoords}
              onChange={e => setPastedCoords(e.target.value)}
              placeholder="Paste SDF / MOL / PDB / XYZ text here..."
              rows={4}
              className="w-full bg-slate-900 border border-slate-700 rounded p-2 text-[11px] text-slate-200 focus:outline-none focus:border-cyan-500 font-mono"
            />
            <div className="flex justify-end">
              <button
                onClick={handlePastedCoordinatesSubmit}
                className="px-3 py-1 bg-cyan-600 hover:bg-cyan-500 text-slate-950 font-bold text-xs rounded transition"
              >
                Parse 3D Coordinates & Cluster
              </button>
            </div>
          </div>
        )}
      </div>

      {/* Main Workbench Layout */}
      <div className="flex-1 flex flex-col lg:flex-row overflow-hidden relative">
        {/* Left / Center: Interactive 3D Visualizer Viewports */}
        <div className="flex-1 flex flex-col min-w-0 bg-slate-950 relative border-r border-slate-800">
          {/* Viewport Control Bar */}
          <div className="h-9 bg-slate-900/90 border-b border-slate-800 px-3 flex items-center justify-between text-xs shrink-0 z-10">
            <div className="flex items-center gap-1.5">
              <button
                onClick={() => setActiveTab3D('both')}
                className={`px-2.5 py-1 rounded text-[11px] transition ${
                  activeTab3D === 'both' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/40 font-bold' : 'text-slate-400 hover:text-white'
                }`}
              >
                Dual View
              </button>
              <button
                onClick={() => setActiveTab3D('pca')}
                className={`px-2.5 py-1 rounded text-[11px] transition ${
                  activeTab3D === 'pca' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/40 font-bold' : 'text-slate-400 hover:text-white'
                }`}
              >
                3D Chemical Space (PCA)
              </button>
              <button
                onClick={() => setActiveTab3D('conformer')}
                className={`px-2.5 py-1 rounded text-[11px] transition ${
                  activeTab3D === 'conformer' ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/40 font-bold' : 'text-slate-400 hover:text-white'
                }`}
              >
                3D Atomic Conformer
              </button>
            </div>

            {/* Clustering K Selector */}
            <div className="flex items-center gap-2">
              <span className="text-[11px] text-slate-400">K-Means:</span>
              <div className="flex items-center gap-1">
                {[2, 3, 4, 5].map(k => (
                  <button
                    key={k}
                    onClick={() => handleKChange(k)}
                    className={`w-5 h-5 rounded text-[10px] font-bold flex items-center justify-center transition ${
                      kClusters === k
                        ? 'bg-cyan-500 text-slate-950 shadow'
                        : 'bg-slate-800 text-slate-400 hover:text-white'
                    }`}
                    title={k === 2 ? 'Silhouette-optimal k=2 (0.479) from paper' : `k=${k}`}
                  >
                    {k}
                  </button>
                ))}
              </div>
              <span className="text-[10px] text-emerald-400 border-l border-slate-800 pl-2">
                Silh: {clusterData.silhouetteScore}
              </span>
            </div>
          </div>

          {/* 3D Visualizer Canvas Grid */}
          <div className="flex-1 relative flex flex-col md:flex-row overflow-hidden">
            {/* 3D Chemical Space PCA */}
            {(activeTab3D === 'both' || activeTab3D === 'pca') && (
              <div className="flex-1 relative h-full min-h-[300px]">
                <PCAViewer3D
                  molecules={library}
                  clusters={clusterData.clusters}
                  activeMolecule={activeRecord}
                  onSelectMolecule={handleSelectMolecule}
                  nearestNeighbors={nearestNeighbors}
                />
              </div>
            )}

            {/* 3D Atomic Conformer Viewer */}
            {(activeTab3D === 'both' || activeTab3D === 'conformer') && (
              <div className="flex-1 relative h-full min-h-[300px]">
                <Conformer3DViewer
                  conformer={activeConformer}
                  compoundName={compoundName}
                />
              </div>
            )}
          </div>
        </div>

        {/* Right Inspector Sidebar: Compound Profile, 7 Descriptors, Nearest Neighbors */}
        <div className="w-full lg:w-96 bg-slate-900 flex flex-col overflow-y-auto p-3 gap-3 shrink-0 border-l border-slate-800">
          {/* Active Compound Badge */}
          <div className="bg-slate-950 p-3 rounded-lg border border-slate-800 space-y-2">
            <div className="flex items-center justify-between">
              <span className="text-[10px] text-cyan-400 uppercase tracking-wider font-semibold">Active Compound</span>
              {activeRecord?.approved && (
                <span className="text-[9px] px-1.5 py-0.2 bg-emerald-500/20 text-emerald-300 border border-emerald-500/40 rounded">
                  FDA Approved
                </span>
              )}
            </div>
            <div className="text-base font-bold text-white truncate">{compoundName}</div>
            <div className="text-[11px] text-slate-400 break-all bg-slate-900 p-1.5 rounded border border-slate-800/80">
              {smilesInput}
            </div>

            <div className="grid grid-cols-2 gap-2 pt-1 text-[11px]">
              <div className="bg-slate-900/60 p-2 rounded border border-slate-800">
                <div className="text-[10px] text-slate-500">Cluster Assignment</div>
                <div className="font-semibold text-cyan-300 flex items-center gap-1.5 mt-0.5">
                  <span
                    className="w-2.5 h-2.5 rounded-full"
                    style={{
                      backgroundColor: clusterData.clusters.find(c => c.id === activeRecord?.cluster)?.color || '#06b6d4'
                    }}
                  />
                  <span>Cluster #{((activeRecord?.cluster || 0) + 1)}</span>
                </div>
              </div>
              <div className="bg-slate-900/60 p-2 rounded border border-slate-800">
                <div className="text-[10px] text-slate-500">3D PCA Coordinates</div>
                <div className="font-semibold text-slate-200 text-[10px] mt-0.5">
                  ({activeRecord?.pca[0] || 0}, {activeRecord?.pca[1] || 0}, {activeRecord?.pca[2] || 0})
                </div>
              </div>
            </div>
          </div>

          {/* 7 RDKit Descriptors Panel */}
          <DescriptorPanel
            descriptors={activeDescriptors}
            drugScore={drugScoreBreakdown}
            compoundName={compoundName}
          />

          {/* Nearest Chemical Space Neighbors */}
          <div className="bg-slate-950 p-3 rounded-lg border border-slate-800 space-y-2 font-mono text-xs">
            <div className="flex items-center justify-between text-[11px] font-semibold text-slate-300">
              <span>Nearest DrugBank Neighbors</span>
              <span className="text-[10px] text-slate-500">Euclidean in R³ PCA</span>
            </div>
            <div className="space-y-1.5">
              {nearestNeighbors.map((nn, i) => (
                <div
                  key={nn.id}
                  onClick={() => {
                    const mol = library.find(m => m.id === nn.id);
                    if (mol) handleSelectMolecule(mol);
                  }}
                  className="p-1.5 bg-slate-900/80 hover:bg-slate-850 rounded border border-slate-800 flex items-center justify-between cursor-pointer transition"
                >
                  <div>
                    <div className="font-semibold text-slate-200 text-[11px] flex items-center gap-1">
                      <span>{nn.name}</span>
                      {nn.approved && (
                        <span className="text-[8px] px-1 bg-emerald-950 text-emerald-400 rounded">FDA</span>
                      )}
                    </div>
                    <div className="text-[10px] text-slate-500">{nn.id}</div>
                  </div>
                  <div className="text-right">
                    <div className="text-[10px] text-cyan-400">Δ = {nn.distance}</div>
                    <div className="text-[9px] text-slate-500">SBOF: {nn.sbofScore.toFixed(3)}</div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>

      {/* Bottom Panel: Top-Ranked Table */}
      <div className="h-64 bg-slate-950 border-t border-slate-800 shrink-0 z-10">
        <TopRankedTable
          molecules={library}
          clusters={clusterData.clusters}
          activeMoleculeId={activeRecord?.id || ''}
          onSelectMolecule={handleSelectMolecule}
        />
      </div>

      {/* Pipeline Integration Modal */}
      <PipelineIntegrationModal
        isOpen={isPipelineModalOpen}
        onClose={() => setIsPipelineModalOpen(false)}
        payload={pipelinePayload}
        conformer={activeConformer}
      />

      {/* 27-Point Dataset Validation Checks Modal */}
      <ValidationChecksModal
        isOpen={isValidationsModalOpen}
        onClose={() => setIsValidationsModalOpen(false)}
        dataset={library}
      />

      {/* Solving Overlay Banner */}
      {isSolving && (
        <div className="fixed bottom-4 left-4 z-50 bg-slate-900/95 border border-cyan-500 text-slate-200 px-4 py-3 rounded-lg shadow-2xl flex items-center gap-3 backdrop-blur font-mono text-xs animate-in slide-in-from-bottom duration-200">
          <RotateCw className="w-4 h-4 text-cyan-400 animate-spin" />
          <div>
            <div className="font-bold text-cyan-300">SBOF Dual-Ascent Engine</div>
            <div className="text-[11px] text-slate-400">{solvingMessage}</div>
          </div>
        </div>
      )}
    </div>
  );
}
