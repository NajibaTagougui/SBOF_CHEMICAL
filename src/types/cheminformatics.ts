/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

export interface MolecularDescriptors {
  mw: number;        // Molecular Weight (Da) [100.0 - 894.1]
  logP: number;      // Octanol-water partition coefficient [-15.0 - 15.0]
  qed: number;       // Quantitative Estimate of Drug-Likeness [0.0 - 1.0]
  hba: number;       // Hydrogen Bond Acceptors (count)
  hbd: number;       // Hydrogen Bond Donors (count)
  psa: number;       // Polar Surface Area (Å²) [0.0 - 542.0]
  rtb: number;       // Rotatable Bonds (count)
}

export interface DrugScoreBreakdown {
  compositeScore: number; // s_i in [0, 1]
  lipinskiPassed: boolean;
  lipinskiViolations: number;
  veberPassed: boolean;
  veberViolations: number;
  sbofScore: number;     // Learned SBOF probability or SBOFNet potential
  sbofRank: number;      // Rank in chemical library
  percentile: number;    // Top X%
}

export interface Atom3D {
  id: number;
  element: string;
  symbol: string;
  x: number;
  y: number;
  z: number;
  formalCharge?: number;
  aromatic?: boolean;
}

export interface Bond3D {
  source: number;
  target: number;
  order: number; // 1: single, 2: double, 3: triple, 1.5: aromatic
}

export interface Conformer3D {
  atoms: Atom3D[];
  bonds: Bond3D[];
  radiusOfGyration: number;
  dimensions: { x: number; y: number; z: number };
  molecularFormula: string;
}

export interface MoleculeRecord {
  id: string;              // DrugBank ID or Compound ID
  name: string;
  smiles: string;
  approved: boolean;       // FDA Approved
  descriptors: MolecularDescriptors;
  pca: [number, number, number]; // PC1, PC2, PC3
  cluster: number;         // 0 to k-1
  sbofScore: number;       // SBOF distribution mass / score
  drugScore: number;       // Composite drug score s_i
  isUserQuery?: boolean;
}

export interface ClusterInfo {
  id: number;
  name: string;
  color: string;
  count: number;
  meanMw: number;
  meanDrugScore: number;
  meanQed: number;
  centroid: [number, number, number];
}

export interface PipelineExportPayload {
  version: string;
  timestamp: string;
  engine: string;
  source: string;
  queryCompound: {
    name: string;
    smiles: string;
    descriptors: MolecularDescriptors;
    pcaCoordinates: { pc1: number; pc2: number; pc3: number };
    clusterAssignment: {
      clusterId: number;
      clusterName: string;
      silhouetteOptimalK: number;
    };
    sbofMetrics: {
      compositeScore: number;
      sbofPotential: number;
      libraryRank: number;
      totalLibrarySize: number;
      percentile: number;
    };
    druglikenessRules: {
      lipinskiRuleOfFive: { passed: boolean; violations: number };
      veberOralBioavailability: { passed: boolean; violations: number };
    };
    nearestNeighbors: Array<{
      id: string;
      name: string;
      distance: number;
      sbofScore: number;
      approved: boolean;
    }>;
  };
}
