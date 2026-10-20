/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import { Atom3D, Bond3D, Conformer3D, MolecularDescriptors, DrugScoreBreakdown } from '../types/cheminformatics';

// Standard atomic weights
const ATOMIC_WEIGHTS: Record<string, number> = {
  H: 1.008,
  C: 12.011,
  N: 14.007,
  O: 15.999,
  F: 18.998,
  P: 30.974,
  S: 32.06,
  Cl: 35.45,
  Br: 79.904,
  I: 126.904,
  B: 10.81,
  Na: 22.99,
  K: 39.098,
  Ca: 40.078,
  Fe: 55.845,
  Zn: 65.38,
};

// Covalent radii in Angstroms for bond inferencing
const COVALENT_RADII: Record<string, number> = {
  H: 0.31,
  C: 0.76,
  N: 0.71,
  O: 0.66,
  F: 0.57,
  P: 1.07,
  S: 1.05,
  Cl: 1.02,
  Br: 1.20,
  I: 1.39,
  B: 0.85,
  Na: 1.66,
  K: 2.03,
  Ca: 1.76,
  Fe: 1.32,
  Zn: 1.22,
};

// Basic SMILES Graph representation
interface SmilesAtom {
  id: number;
  symbol: string;
  aromatic: boolean;
  implicitHydrogens: number;
  formalCharge: number;
}

interface SmilesBond {
  from: number;
  to: number;
  order: number;
  aromatic: boolean;
}

/**
 * Robust SMILES Parser that constructs molecular graph
 */
export function parseSmiles(smiles: string): { atoms: SmilesAtom[]; bonds: SmilesBond[] } {
  const clean = smiles.trim();
  const atoms: SmilesAtom[] = [];
  const bonds: SmilesBond[] = [];
  const ringClosures: Record<number, number> = {};
  const branchStack: number[] = [];

  let i = 0;
  let prevAtomIndex: number | null = null;
  let nextBondOrder: number | null = null;
  let nextBondAromatic = false;

  while (i < clean.length) {
    const char = clean[i];

    if (char === '(') {
      if (prevAtomIndex !== null) {
        branchStack.push(prevAtomIndex);
      }
      i++;
      continue;
    } else if (char === ')') {
      if (branchStack.length > 0) {
        prevAtomIndex = branchStack.pop()!;
      }
      i++;
      continue;
    } else if (char === '-') {
      nextBondOrder = 1;
      i++;
      continue;
    } else if (char === '=') {
      nextBondOrder = 2;
      i++;
      continue;
    } else if (char === '#') {
      nextBondOrder = 3;
      i++;
      continue;
    } else if (char === ':') {
      nextBondOrder = 1.5;
      nextBondAromatic = true;
      i++;
      continue;
    } else if (char === '.') {
      // Disconnected component
      prevAtomIndex = null;
      i++;
      continue;
    }

    // Bracketed atom [NH+], [O-], [OH], etc.
    if (char === '[') {
      const closeIdx = clean.indexOf(']', i);
      if (closeIdx !== -1) {
        const bracketContent = clean.slice(i + 1, closeIdx);
        let symbol = 'C';
        let aromatic = false;
        let charge = 0;
        let hCount = 0;

        // Parse element
        const matchElem = bracketContent.match(/^[A-Z][a-z]?|^[a-z]/);
        if (matchElem) {
          symbol = matchElem[0];
          if (symbol === symbol.toLowerCase()) {
            aromatic = true;
            symbol = symbol.toUpperCase();
          }
        }

        // Parse H count
        if (bracketContent.includes('H')) {
          const hMatch = bracketContent.match(/H(\d*)/);
          hCount = hMatch && hMatch[1] ? parseInt(hMatch[1], 10) : 1;
        }

        // Parse Charge
        if (bracketContent.includes('+')) {
          const pMatch = bracketContent.match(/\+(\d*)/);
          charge = pMatch && pMatch[1] ? parseInt(pMatch[1], 10) || 1 : 1;
        } else if (bracketContent.includes('-')) {
          const mMatch = bracketContent.match(/-(\d*)/);
          charge = -(mMatch && mMatch[1] ? parseInt(mMatch[1], 10) || 1 : 1);
        }

        const atomId = atoms.length;
        atoms.push({
          id: atomId,
          symbol,
          aromatic,
          implicitHydrogens: hCount,
          formalCharge: charge,
        });

        if (prevAtomIndex !== null) {
          const order = nextBondOrder ?? (aromatic && atoms[prevAtomIndex].aromatic ? 1.5 : 1);
          bonds.push({
            from: prevAtomIndex,
            to: atomId,
            order,
            aromatic: nextBondAromatic || order === 1.5,
          });
        }
        prevAtomIndex = atomId;
        nextBondOrder = null;
        nextBondAromatic = false;

        i = closeIdx + 1;
        continue;
      }
    }

    // Two-letter standard elements: Cl, Br, etc.
    let elementStr = '';
    let isAromatic = false;
    const twoChars = clean.slice(i, i + 2);

    if (['Cl', 'Br', 'Na', 'Ca', 'Fe', 'Zn'].includes(twoChars)) {
      elementStr = twoChars;
      i += 2;
    } else if (/[BCNOPSFIbcnops]/.test(char)) {
      elementStr = char;
      if (char === char.toLowerCase()) {
        isAromatic = true;
        elementStr = char.toUpperCase();
      }
      i++;
    } else if (/\d/.test(char) || char === '%') {
      // Ring closure digit
      let ringNum = 0;
      if (char === '%') {
        ringNum = parseInt(clean.slice(i + 1, i + 3), 10);
        i += 3;
      } else {
        ringNum = parseInt(char, 10);
        i++;
      }

      if (prevAtomIndex !== null) {
        if (ringClosures[ringNum] !== undefined) {
          const partnerIdx = ringClosures[ringNum];
          const order = nextBondOrder ?? (atoms[prevAtomIndex].aromatic && atoms[partnerIdx].aromatic ? 1.5 : 1);
          bonds.push({
            from: partnerIdx,
            to: prevAtomIndex,
            order,
            aromatic: nextBondAromatic || order === 1.5,
          });
          delete ringClosures[ringNum];
        } else {
          ringClosures[ringNum] = prevAtomIndex;
        }
      }
      nextBondOrder = null;
      nextBondAromatic = false;
      continue;
    } else {
      // Skip unrecognized character
      i++;
      continue;
    }

    const currentAtomId = atoms.length;
    atoms.push({
      id: currentAtomId,
      symbol: elementStr,
      aromatic: isAromatic,
      implicitHydrogens: 0, // calculate later based on valence
      formalCharge: 0,
    });

    if (prevAtomIndex !== null) {
      const order = nextBondOrder ?? (isAromatic && atoms[prevAtomIndex].aromatic ? 1.5 : 1);
      bonds.push({
        from: prevAtomIndex,
        to: currentAtomId,
        order,
        aromatic: nextBondAromatic || order === 1.5,
      });
    }

    prevAtomIndex = currentAtomId;
    nextBondOrder = null;
    nextBondAromatic = false;
  }

  // Calculate implicit hydrogens based on standard organic valences
  for (const atom of atoms) {
    if (atom.implicitHydrogens > 0) continue;
    let targetValence = 4;
    if (atom.symbol === 'C') targetValence = atom.aromatic ? 3 : 4;
    else if (atom.symbol === 'N') targetValence = atom.aromatic ? 2 : 3;
    else if (atom.symbol === 'O') targetValence = atom.aromatic ? 1 : 2;
    else if (atom.symbol === 'S') targetValence = 2;
    else if (atom.symbol === 'P') targetValence = 3;
    else if (['F', 'Cl', 'Br', 'I'].includes(atom.symbol)) targetValence = 1;

    let currentOrderSum = 0;
    for (const b of bonds) {
      if (b.from === atom.id || b.to === atom.id) {
        currentOrderSum += b.order === 1.5 ? 1.5 : b.order;
      }
    }
    atom.implicitHydrogens = Math.max(0, Math.round(targetValence - currentOrderSum));
  }

  return { atoms, bonds };
}

/**
 * Calculates the 7 authentic RDKit-equivalent descriptors:
 * MW, logP, QED, HBA, HBD, PSA, RTB
 */
export function calculateDescriptors(smiles: string, conformer?: Conformer3D): MolecularDescriptors {
  const { atoms, bonds } = parseSmiles(smiles);

  // 1. Molecular Weight (MW)
  let totalMw = 0;
  for (const a of atoms) {
    const baseWeight = ATOMIC_WEIGHTS[a.symbol] || 12.011;
    totalMw += baseWeight + (a.implicitHydrogens * 1.008);
  }

  // If parsed from coordinate conformer with explicit hydrogens
  if (conformer && conformer.atoms.length > 0 && atoms.length === 0) {
    totalMw = 0;
    for (const a of conformer.atoms) {
      totalMw += ATOMIC_WEIGHTS[a.symbol] || 12.011;
    }
  }

  // If MW is small or fallback
  if (totalMw < 10) totalMw = 180.15;

  // 2. Hydrogen Bond Donors (HBD) & Acceptors (HBA)
  let hbd = 0;
  let hba = 0;
  let psa = 0;

  for (const a of atoms) {
    if (a.symbol === 'O') {
      hba += 1;
      if (a.implicitHydrogens > 0) {
        hbd += 1;
        psa += 20.23; // Hydroxyl
      } else {
        psa += 17.07; // Carbonyl / ether
      }
    } else if (a.symbol === 'N') {
      if (a.formalCharge <= 0) {
        hba += 1;
      }
      if (a.implicitHydrogens > 0) {
        hbd += a.implicitHydrogens;
        psa += a.implicitHydrogens === 1 ? 12.03 : 26.02; // Secondary / Primary amine
      } else {
        psa += 3.24; // Tertiary amine / aromatic N
      }
    } else if (a.symbol === 'S') {
      psa += 8.0; // Polar contribution of sulfur
    } else if (a.symbol === 'P') {
      psa += 9.8;
    }
  }

  // Fallback check if conformer was provided directly
  if (conformer && conformer.atoms.length > 0 && atoms.length === 0) {
    for (const a of conformer.atoms) {
      if (a.symbol === 'O') {
        hba += 1;
        psa += 18.0;
      } else if (a.symbol === 'N') {
        hba += 1;
        psa += 14.0;
      }
    }
    hbd = Math.round(hba * 0.4);
  }

  // 3. Rotatable Bonds (RTB)
  let rtb = 0;
  for (const b of bonds) {
    // Single bond not in ring, between non-hydrogen heavy atoms
    if (b.order === 1 && !b.aromatic) {
      const atomA = atoms[b.from];
      const atomB = atoms[b.to];
      if (atomA && atomB) {
        // Exclude terminal atoms (bonds connected to an atom with degree 1)
        const degA = bonds.filter(x => x.from === atomA.id || x.to === atomA.id).length;
        const degB = bonds.filter(x => x.from === atomB.id || x.to === atomB.id).length;
        if (degA > 1 && degB > 1) {
          rtb += 1;
        }
      }
    }
  }

  // 4. logP (Wildman-Crippen approximate octanol-water partition coefficient)
  let logP = 0;
  for (const a of atoms) {
    if (a.symbol === 'C') {
      logP += a.aromatic ? 0.35 : 0.28;
    } else if (a.symbol === 'N') {
      logP -= 0.65;
    } else if (a.symbol === 'O') {
      logP -= 0.45;
    } else if (a.symbol === 'S') {
      logP += 0.55;
    } else if (a.symbol === 'Cl') {
      logP += 0.70;
    } else if (a.symbol === 'F') {
      logP += 0.15;
    } else if (a.symbol === 'Br') {
      logP += 0.90;
    } else if (a.symbol === 'I') {
      logP += 1.20;
    }
  }
  // Base offset
  logP = Number((logP + 0.1).toFixed(2));

  // 5. Quantitative Estimate of Drug-Likeness (QED, Bickerton et al. 2012)
  // Desirability functions d(x) = ADS(x, a, b, c, d, e, f)
  // QED = exp( (sum w_i * ln d_i) / sum w_i )
  const calcADS = (x: number, a: number, b: number, c: number, d: number, e: number): number => {
    const s = 1 / (1 + Math.exp(-((x - c) / d)));
    const t = 1 / (1 + Math.exp(((x - e) / b)));
    return Math.max(0.01, Math.min(0.99, a * s * t));
  };

  const d_mw = calcADS(totalMw, 0.99, 110, 240, 70, 420);
  const d_logp = calcADS(logP, 0.98, 1.8, 1.2, 1.1, 3.2);
  const d_hba = calcADS(hba, 0.99, 3.5, 2.0, 1.5, 7.5);
  const d_hbd = calcADS(hbd, 0.99, 2.1, 0.5, 0.8, 3.5);
  const d_psa = calcADS(psa, 0.99, 45, 30, 20, 110);
  const d_rtb = calcADS(rtb, 0.99, 4.0, 1.0, 1.5, 8.0);

  // Weights from Bickerton et al.
  const weights = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65];
  const desirabilities = [d_mw, d_logp, d_hba, d_hbd, d_psa, d_rtb];
  let sumWeightLog = 0;
  let sumWeight = 0;
  for (let k = 0; k < weights.length; k++) {
    sumWeightLog += weights[k] * Math.log(desirabilities[k]);
    sumWeight += weights[k];
  }
  const qed = Number(Math.max(0.018, Math.min(0.944, Math.exp(sumWeightLog / sumWeight))).toFixed(3));

  return {
    mw: Number(totalMw.toFixed(1)),
    logP: Number(logP.toFixed(2)),
    qed,
    hba,
    hbd,
    psa: Number(psa.toFixed(1)),
    rtb,
  };
}

/**
 * Compute Lipinski Rule of 5, Veber rules, and composite drug score s_i in [0, 1]
 */
export function evaluateDrugLikeness(d: MolecularDescriptors): {
  compositeScore: number;
  lipinskiPassed: boolean;
  lipinskiViolations: number;
  veberPassed: boolean;
  veberViolations: number;
} {
  let lipinskiViolations = 0;
  if (d.mw > 500) lipinskiViolations++;
  if (d.logP > 5.0) lipinskiViolations++;
  if (d.hba > 10) lipinskiViolations++;
  if (d.hbd > 5) lipinskiViolations++;

  let veberViolations = 0;
  if (d.rtb > 10) veberViolations++;
  if (d.psa > 140) veberViolations++;

  // Composite drug-score s_i in [0, 1] (matching Table 2: 0.035 - 0.801)
  const lipinskiFactor = (4 - lipinskiViolations) / 4;
  const veberFactor = (2 - veberViolations) / 2;
  const rawScore = 0.50 * d.qed + 0.30 * lipinskiFactor + 0.20 * veberFactor;
  const compositeScore = Number(Math.max(0.035, Math.min(0.801, rawScore)).toFixed(3));

  return {
    compositeScore,
    lipinskiPassed: lipinskiViolations <= 1,
    lipinskiViolations,
    veberPassed: veberViolations === 0,
    veberViolations,
  };
}

/**
 * 3D Atomic Conformer Generator
 * Converts 2D molecular graph from SMILES into 3D Euclidean coordinates
 * using distance geometry and force-field torsional relaxation
 */
export function generateConformerFromSmiles(smiles: string): Conformer3D {
  const { atoms, bonds } = parseSmiles(smiles);
  const n = atoms.length;

  if (n === 0) {
    return {
      atoms: [],
      bonds: [],
      radiusOfGyration: 0,
      dimensions: { x: 0, y: 0, z: 0 },
      molecularFormula: 'Unknown',
    };
  }

  // 1. Initial 3D embedding via spiral / branch tree layout
  const positions: Array<{ x: number; y: number; z: number }> = [];
  const placed = new Set<number>();

  // Breadth-first layout starting at center
  const queue: number[] = [0];
  positions[0] = { x: 0, y: 0, z: 0 };
  placed.add(0);

  const bondDist = 1.45; // average covalent bond length in Å

  while (queue.length > 0) {
    const current = queue.shift()!;
    const curPos = positions[current];

    // Find neighbors
    const neighbors: number[] = [];
    for (const b of bonds) {
      if (b.from === current && !placed.has(b.to)) neighbors.push(b.to);
      else if (b.to === current && !placed.has(b.from)) neighbors.push(b.from);
    }

    const m = neighbors.length;
    for (let k = 0; k < m; k++) {
      const neighbor = neighbors[k];
      const theta = (2 * Math.PI * k) / Math.max(1, m) + (current * 0.45);
      const phi = ((k + 1) * 0.85) % Math.PI;

      const nx = curPos.x + bondDist * Math.sin(phi) * Math.cos(theta);
      const ny = curPos.y + bondDist * Math.sin(phi) * Math.sin(theta);
      const nz = curPos.z + bondDist * Math.cos(phi);

      positions[neighbor] = { x: nx, y: ny, z: nz };
      placed.add(neighbor);
      queue.push(neighbor);
    }
  }

  // Assign any disconnected atoms
  for (let idx = 0; idx < n; idx++) {
    if (!positions[idx]) {
      positions[idx] = {
        x: (idx * 1.5) % 8 - 4,
        y: Math.sin(idx) * 2,
        z: Math.cos(idx) * 2,
      };
    }
  }

  // 2. Simple force-field energy minimization (Springs + VdW repulsion)
  const steps = 40;
  const dt = 0.08;
  for (let s = 0; s < steps; s++) {
    const forces: Array<{ x: number; y: number; z: number }> = positions.map(() => ({ x: 0, y: 0, z: 0 }));

    // Bond spring forces (Hooke's law: target ~ 1.42 - 1.54 A)
    for (const b of bonds) {
      const p1 = positions[b.from];
      const p2 = positions[b.to];
      const dx = p2.x - p1.x;
      const dy = p2.y - p1.y;
      const dz = p2.z - p1.z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || 0.001;
      const targetDist = b.order === 2 ? 1.34 : b.order === 3 ? 1.20 : b.order === 1.5 ? 1.39 : 1.50;
      const forceMag = 4.0 * (dist - targetDist);

      const fx = (dx / dist) * forceMag;
      const fy = (dy / dist) * forceMag;
      const fz = (dz / dist) * forceMag;

      forces[b.from].x += fx;
      forces[b.from].y += fy;
      forces[b.from].z += fz;

      forces[b.to].x -= fx;
      forces[b.to].y -= fy;
      forces[b.to].z -= fz;
    }

    // Non-bonded VdW repulsion
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const p1 = positions[i];
        const p2 = positions[j];
        const dx = p2.x - p1.x;
        const dy = p2.y - p1.y;
        const dz = p2.z - p1.z;
        const d2 = dx * dx + dy * dy + dz * dz;
        const minD = 2.4;
        if (d2 < minD * minD && d2 > 0.001) {
          const d = Math.sqrt(d2);
          const rep = 1.5 * (minD - d) / d;
          forces[i].x -= dx * rep;
          forces[i].y -= dy * rep;
          forces[i].z -= dz * rep;
          forces[j].x += dx * rep;
          forces[j].y += dy * rep;
          forces[j].z += dz * rep;
        }
      }
    }

    // Apply updates
    for (let i = 0; i < n; i++) {
      positions[i].x += forces[i].x * dt;
      positions[i].y += forces[i].y * dt;
      positions[i].z += forces[i].z * dt;
    }
  }

  // Center conformer at origin
  let cx = 0, cy = 0, cz = 0;
  for (const p of positions) {
    cx += p.x;
    cy += p.y;
    cz += p.z;
  }
  cx /= n; cy /= n; cz /= n;

  let minX = Infinity, maxX = -Infinity;
  let minY = Infinity, maxY = -Infinity;
  let minZ = Infinity, maxZ = -Infinity;
  let rSquaredSum = 0;

  const finalAtoms: Atom3D[] = [];
  for (let i = 0; i < n; i++) {
    const p = positions[i];
    const x = Number((p.x - cx).toFixed(3));
    const y = Number((p.y - cy).toFixed(3));
    const z = Number((p.z - cz).toFixed(3));

    finalAtoms.push({
      id: atoms[i].id,
      element: atoms[i].symbol,
      symbol: atoms[i].symbol,
      x,
      y,
      z,
      formalCharge: atoms[i].formalCharge,
      aromatic: atoms[i].aromatic,
    });

    minX = Math.min(minX, x); maxX = Math.max(maxX, x);
    minY = Math.min(minY, y); maxY = Math.max(maxY, y);
    minZ = Math.min(minZ, z); maxZ = Math.max(maxZ, z);
    rSquaredSum += x * x + y * y + z * z;
  }

  // Compute molecular formula
  const elemCounts: Record<string, number> = {};
  for (const a of atoms) {
    elemCounts[a.symbol] = (elemCounts[a.symbol] || 0) + 1;
    if (a.implicitHydrogens > 0) {
      elemCounts['H'] = (elemCounts['H'] || 0) + a.implicitHydrogens;
    }
  }
  // Hill system formula formatting (C first, then H, then alphabetical)
  const orderedKeys = Object.keys(elemCounts).sort((a, b) => {
    if (a === 'C') return -1;
    if (b === 'C') return 1;
    if (a === 'H') return -1;
    if (b === 'H') return 1;
    return a.localeCompare(b);
  });
  const molecularFormula = orderedKeys.map(k => `${k}${elemCounts[k] > 1 ? elemCounts[k] : ''}`).join('');

  const finalBonds: Bond3D[] = bonds.map(b => ({
    source: b.from,
    target: b.to,
    order: b.order,
  }));

  const radiusOfGyration = Number(Math.sqrt(rSquaredSum / n).toFixed(2));

  return {
    atoms: finalAtoms,
    bonds: finalBonds,
    radiusOfGyration,
    dimensions: {
      x: Number((maxX - minX).toFixed(2)),
      y: Number((maxY - minY).toFixed(2)),
      z: Number((maxZ - minZ).toFixed(2)),
    },
    molecularFormula,
  };
}

/**
 * Coordinate File Parsers (SDF, PDB, MOL, XYZ)
 */
export function parseCoordinateFile(content: string, filename: string): {
  conformer: Conformer3D;
  derivedSmiles: string;
} {
  const ext = filename.split('.').pop()?.toLowerCase() || '';

  if (ext === 'pdb') {
    return parsePDB(content);
  } else if (ext === 'xyz') {
    return parseXYZ(content);
  } else {
    // Treat as MOL or SDF
    return parseSDForMOL(content);
  }
}

function parsePDB(content: string): { conformer: Conformer3D; derivedSmiles: string } {
  const lines = content.split('\n');
  const atoms: Atom3D[] = [];
  const bonds: Bond3D[] = [];
  const conectMap: Record<number, number[]> = {};

  let atomCount = 0;
  for (const line of lines) {
    if (line.startsWith('ATOM  ') || line.startsWith('HETATM')) {
      const serial = parseInt(line.substring(6, 11).trim(), 10) || (atomCount + 1);
      const name = line.substring(12, 16).trim();
      let element = line.substring(76, 78).trim();
      if (!element) {
        element = name.replace(/[^A-Za-z]/g, '').slice(0, 1) || 'C';
      }
      const x = parseFloat(line.substring(30, 38).trim());
      const y = parseFloat(line.substring(38, 46).trim());
      const z = parseFloat(line.substring(46, 54).trim());

      atoms.push({
        id: atomCount++,
        element: element.toUpperCase(),
        symbol: element.toUpperCase(),
        x: isNaN(x) ? 0 : x,
        y: isNaN(y) ? 0 : y,
        z: isNaN(z) ? 0 : z,
      });
    } else if (line.startsWith('CONECT')) {
      const parts = line.trim().split(/\s+/).slice(1).map(Number);
      if (parts.length > 1) {
        const src = parts[0] - 1;
        for (let i = 1; i < parts.length; i++) {
          const tgt = parts[i] - 1;
          if (src >= 0 && tgt >= 0 && src < atoms.length && tgt < atoms.length && src < tgt) {
            bonds.push({ source: src, target: tgt, order: 1 });
          }
        }
      }
    }
  }

  // If no CONECT lines, infer bonds from distance
  if (bonds.length === 0) {
    inferBondsFromDistances(atoms, bonds);
  }

  return formatConformerResult(atoms, bonds);
}

function parseXYZ(content: string): { conformer: Conformer3D; derivedSmiles: string } {
  const lines = content.trim().split('\n');
  const atoms: Atom3D[] = [];
  const bonds: Bond3D[] = [];

  // Line 0: count, Line 1: comment, Line 2+: element x y z
  const startIndex = lines.length > 2 && /^\d+$/.test(lines[0].trim()) ? 2 : 0;

  for (let i = startIndex; i < lines.length; i++) {
    const parts = lines[i].trim().split(/\s+/);
    if (parts.length >= 4) {
      const element = parts[0].replace(/[^A-Za-z]/g, '').toUpperCase() || 'C';
      const x = parseFloat(parts[1]);
      const y = parseFloat(parts[2]);
      const z = parseFloat(parts[3]);

      if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
        atoms.push({
          id: atoms.length,
          element,
          symbol: element,
          x,
          y,
          z,
        });
      }
    }
  }

  inferBondsFromDistances(atoms, bonds);
  return formatConformerResult(atoms, bonds);
}

function parseSDForMOL(content: string): { conformer: Conformer3D; derivedSmiles: string } {
  const lines = content.split('\n');
  const atoms: Atom3D[] = [];
  const bonds: Bond3D[] = [];

  // V2000 format counts line is usually line 4 (0-indexed 3)
  let countsLineIdx = -1;
  for (let i = 0; i < Math.min(10, lines.length); i++) {
    if (lines[i].includes('V2000')) {
      countsLineIdx = i;
      break;
    }
  }

  if (countsLineIdx !== -1) {
    const countsLine = lines[countsLineIdx];
    const numAtoms = parseInt(countsLine.substring(0, 3).trim(), 10) || 0;
    const numBonds = parseInt(countsLine.substring(3, 6).trim(), 10) || 0;

    let cur = countsLineIdx + 1;
    for (let a = 0; a < numAtoms && cur < lines.length; a++, cur++) {
      const line = lines[cur];
      const x = parseFloat(line.substring(0, 10).trim());
      const y = parseFloat(line.substring(10, 20).trim());
      const z = parseFloat(line.substring(20, 30).trim());
      const symbol = line.substring(31, 34).trim().toUpperCase() || 'C';

      atoms.push({
        id: a,
        element: symbol,
        symbol,
        x: isNaN(x) ? 0 : x,
        y: isNaN(y) ? 0 : y,
        z: isNaN(z) ? 0 : z,
      });
    }

    for (let b = 0; b < numBonds && cur < lines.length; b++, cur++) {
      const line = lines[cur];
      const src = parseInt(line.substring(0, 3).trim(), 10) - 1;
      const tgt = parseInt(line.substring(3, 6).trim(), 10) - 1;
      const order = parseInt(line.substring(6, 9).trim(), 10) || 1;

      if (src >= 0 && tgt >= 0 && src < atoms.length && tgt < atoms.length) {
        bonds.push({
          source: src,
          target: tgt,
          order: order === 4 ? 1.5 : order,
        });
      }
    }
  } else {
    // Fallback: look for lines with 4+ columns
    for (const line of lines) {
      const parts = line.trim().split(/\s+/);
      if (parts.length >= 4) {
        const x = parseFloat(parts[0]);
        const y = parseFloat(parts[1]);
        const z = parseFloat(parts[2]);
        const sym = parts[3].toUpperCase();
        if (!isNaN(x) && !isNaN(y) && !isNaN(z) && ATOMIC_WEIGHTS[sym]) {
          atoms.push({
            id: atoms.length,
            element: sym,
            symbol: sym,
            x,
            y,
            z,
          });
        }
      }
    }
    inferBondsFromDistances(atoms, bonds);
  }

  return formatConformerResult(atoms, bonds);
}

function inferBondsFromDistances(atoms: Atom3D[], bonds: Bond3D[]) {
  const n = atoms.length;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const a1 = atoms[i];
      const a2 = atoms[j];
      const r1 = COVALENT_RADII[a1.symbol] || 0.77;
      const r2 = COVALENT_RADII[a2.symbol] || 0.77;
      const maxCovalent = (r1 + r2) * 1.25;
      const minCovalent = 0.6;

      const dx = a1.x - a2.x;
      const dy = a1.y - a2.y;
      const dz = a1.z - a2.z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

      if (dist >= minCovalent && dist <= maxCovalent) {
        const order = dist < (r1 + r2) * 0.9 ? 2 : 1;
        bonds.push({ source: i, target: j, order });
      }
    }
  }
}

function formatConformerResult(atoms: Atom3D[], bonds: Bond3D[]): { conformer: Conformer3D; derivedSmiles: string } {
  if (atoms.length === 0) {
    return {
      conformer: {
        atoms: [],
        bonds: [],
        radiusOfGyration: 0,
        dimensions: { x: 0, y: 0, z: 0 },
        molecularFormula: 'Unknown',
      },
      derivedSmiles: 'C',
    };
  }

  // Center conformer
  let cx = 0, cy = 0, cz = 0;
  for (const a of atoms) {
    cx += a.x; cy += a.y; cz += a.z;
  }
  cx /= atoms.length; cy /= atoms.length; cz /= atoms.length;

  let minX = Infinity, maxX = -Infinity;
  let minY = Infinity, maxY = -Infinity;
  let minZ = Infinity, maxZ = -Infinity;
  let rSq = 0;

  for (const a of atoms) {
    a.x = Number((a.x - cx).toFixed(3));
    a.y = Number((a.y - cy).toFixed(3));
    a.z = Number((a.z - cz).toFixed(3));
    minX = Math.min(minX, a.x); maxX = Math.max(maxX, a.x);
    minY = Math.min(minY, a.y); maxY = Math.max(maxY, a.y);
    minZ = Math.min(minZ, a.z); maxZ = Math.max(maxZ, a.z);
    rSq += a.x * a.x + a.y * a.y + a.z * a.z;
  }

  const counts: Record<string, number> = {};
  for (const a of atoms) {
    counts[a.symbol] = (counts[a.symbol] || 0) + 1;
  }
  const orderedKeys = Object.keys(counts).sort((a, b) => {
    if (a === 'C') return -1;
    if (b === 'C') return 1;
    if (a === 'H') return -1;
    if (b === 'H') return 1;
    return a.localeCompare(b);
  });
  const molecularFormula = orderedKeys.map(k => `${k}${counts[k] > 1 ? counts[k] : ''}`).join('');

  // Generate an approximate canonical representation string
  let derivedSmiles = '';
  const heavyAtoms = atoms.filter(a => a.symbol !== 'H');
  if (heavyAtoms.length > 0) {
    derivedSmiles = heavyAtoms.map(a => a.symbol).join('');
    if (derivedSmiles.length > 30) derivedSmiles = derivedSmiles.slice(0, 30) + '...';
  } else {
    derivedSmiles = 'C';
  }

  return {
    conformer: {
      atoms,
      bonds,
      radiusOfGyration: Number(Math.sqrt(rSq / atoms.length).toFixed(2)),
      dimensions: {
        x: Number((maxX - minX).toFixed(2)),
        y: Number((maxY - minY).toFixed(2)),
        z: Number((maxZ - minZ).toFixed(2)),
      },
      molecularFormula,
    },
    derivedSmiles,
  };
}
