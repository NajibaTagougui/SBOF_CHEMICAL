/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { Conformer3D, Atom3D } from '../types/cheminformatics';
import { RotateCcw, Compass, Box, Layers, Maximize2, Minimize2 } from 'lucide-react';

interface Conformer3DViewerProps {
  conformer: Conformer3D | null;
  compoundName?: string;
}

// CPK Element Colors
const ELEMENT_COLORS: Record<string, number> = {
  C: 0x475569, // Slate Carbon
  H: 0xf8fafc, // Pure White
  O: 0xef4444, // Red Oxygen
  N: 0x3b82f6, // Blue Nitrogen
  S: 0xeab308, // Yellow Sulfur
  P: 0xf97316, // Orange Phosphorus
  F: 0x10b981, // Emerald Halogen
  CL: 0x10b981,
  BR: 0xa855f7, // Purple Bromine
  I: 0x7c3aed, // Violet Iodine
  B: 0xf472b6,
  FE: 0xe11d48,
  ZN: 0x64748b,
};

// Approximate Van der Waals radii in Ångströms
const VDW_RADII: Record<string, number> = {
  H: 1.10,
  C: 1.70,
  N: 1.55,
  O: 1.52,
  F: 1.47,
  P: 1.80,
  S: 1.80,
  CL: 1.75,
  BR: 1.85,
  I: 1.98,
};

type RenderStyle = 'ballStick' | 'spacefill' | 'sticks' | 'wireframe';

export const Conformer3DViewer: React.FC<Conformer3DViewerProps> = ({ conformer, compoundName }) => {
  const mountRef = useRef<HTMLDivElement>(null);
  const [style, setStyle] = useState<RenderStyle>('ballStick');
  const [showHydrogens, setShowHydrogens] = useState(true);
  const [autoRotate, setAutoRotate] = useState(false);
  const [hoveredAtom, setHoveredAtom] = useState<Atom3D | null>(null);
  const [tooltipPos, setTooltipPos] = useState<{ x: number; y: number } | null>(null);

  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const moleculeGroupRef = useRef<THREE.Group | null>(null);

  const controlsRef = useRef<{
    isDragging: boolean;
    prevX: number;
    prevY: number;
    rotX: number;
    rotY: number;
    distance: number;
  }>({
    isDragging: false,
    prevX: 0,
    prevY: 0,
    rotX: 0.2,
    rotY: 0.3,
    distance: 12,
  });

  // Setup Three.js scene
  useEffect(() => {
    if (!mountRef.current) return;
    const width = mountRef.current.clientWidth;
    const height = mountRef.current.clientHeight;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x070c18); // Dark scientific viewport
    sceneRef.current = scene;

    const camera = new THREE.PerspectiveCamera(40, width / height, 0.1, 100);
    cameraRef.current = camera;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.shadowMap.enabled = true;
    mountRef.current.innerHTML = '';
    mountRef.current.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Lights
    const ambLight = new THREE.AmbientLight(0xffffff, 0.9);
    scene.add(ambLight);

    const dirLight = new THREE.DirectionalLight(0xffffff, 1.4);
    dirLight.position.set(10, 15, 20);
    scene.add(dirLight);

    const rimLight = new THREE.DirectionalLight(0x38bdf8, 0.7);
    rimLight.position.set(-15, -10, -10);
    scene.add(rimLight);

    const moleculeGroup = new THREE.Group();
    scene.add(moleculeGroup);
    moleculeGroupRef.current = moleculeGroup;

    // Animation loop
    let animId: number;
    const animate = () => {
      animId = requestAnimationFrame(animate);

      if (autoRotate && moleculeGroupRef.current) {
        controlsRef.current.rotY += 0.008;
      }

      const ctrl = controlsRef.current;
      const x = ctrl.distance * Math.sin(ctrl.rotY) * Math.cos(ctrl.rotX);
      const y = ctrl.distance * Math.sin(ctrl.rotX);
      const z = ctrl.distance * Math.cos(ctrl.rotY) * Math.cos(ctrl.rotX);

      camera.position.set(x, y, z);
      camera.lookAt(0, 0, 0);

      renderer.render(scene, camera);
    };
    animate();

    const handleResize = () => {
      if (!mountRef.current || !rendererRef.current || !cameraRef.current) return;
      const w = mountRef.current.clientWidth;
      const h = mountRef.current.clientHeight;
      cameraRef.current.aspect = w / h;
      cameraRef.current.updateProjectionMatrix();
      rendererRef.current.setSize(w, h);
    };
    window.addEventListener('resize', handleResize);

    return () => {
      cancelAnimationFrame(animId);
      window.removeEventListener('resize', handleResize);
      renderer.dispose();
      if (mountRef.current) mountRef.current.innerHTML = '';
    };
  }, [autoRotate]);

  // Build 3D Conformer Geometry
  useEffect(() => {
    if (!moleculeGroupRef.current || !conformer) return;
    const group = moleculeGroupRef.current;

    // Clean previous children
    while (group.children.length > 0) {
      const child = group.children[0];
      group.remove(child);
      if (child instanceof THREE.Mesh || child instanceof THREE.LineSegments) {
        child.geometry.dispose();
        if (Array.isArray(child.material)) child.material.forEach(m => m.dispose());
        else child.material.dispose();
      }
    }

    const atoms = conformer.atoms.filter(a => showHydrogens || a.element !== 'H');
    const atomIndexMap = new Map<number, number>();
    atoms.forEach((a, i) => atomIndexMap.set(a.id, i));

    // Determine scale & bounding size
    const maxDim = Math.max(conformer.dimensions.x, conformer.dimensions.y, conformer.dimensions.z, 2.0);
    controlsRef.current.distance = Math.max(7, maxDim * 2.3);

    // 1. Render Atoms
    atoms.forEach(atom => {
      const sym = atom.symbol.toUpperCase();
      const color = ELEMENT_COLORS[sym] || 0x64748b;
      const vdw = VDW_RADII[sym] || 1.70;

      let radius = 0.38;
      if (style === 'spacefill') {
        radius = vdw * 0.75;
      } else if (style === 'sticks' || style === 'wireframe') {
        radius = 0.08;
      } else {
        // Ball and stick
        radius = atom.element === 'H' ? 0.22 : 0.38;
      }

      const geo = new THREE.SphereGeometry(radius, 20, 20);
      const mat = new THREE.MeshStandardMaterial({
        color,
        roughness: 0.25,
        metalness: 0.2,
      });

      const mesh = new THREE.Mesh(geo, mat);
      mesh.position.set(atom.x, atom.y, atom.z);
      mesh.userData = { atom };
      group.add(mesh);
    });

    // 2. Render Bonds
    if (style !== 'spacefill') {
      conformer.bonds.forEach(bond => {
        const atom1 = atoms.find(a => a.id === bond.source);
        const atom2 = atoms.find(a => a.id === bond.target);
        if (!atom1 || !atom2) return;

        const p1 = new THREE.Vector3(atom1.x, atom1.y, atom1.z);
        const p2 = new THREE.Vector3(atom2.x, atom2.y, atom2.z);
        const dist = p1.distanceTo(p2);
        const midPoint = new THREE.Vector3().addVectors(p1, p2).multiplyScalar(0.5);

        if (style === 'wireframe') {
          const lineGeo = new THREE.BufferGeometry().setFromPoints([p1, p2]);
          const lineMat = new THREE.LineBasicMaterial({ color: 0x94a3b8, linewidth: 2 });
          const line = new THREE.Line(lineGeo, lineMat);
          group.add(line);
        } else {
          // Cylinder bond
          const radius = style === 'sticks' ? 0.14 : 0.10;
          const cylGeo = new THREE.CylinderGeometry(radius, radius, dist, 12);
          const cylMat = new THREE.MeshStandardMaterial({
            color: 0x94a3b8,
            roughness: 0.4,
            metalness: 0.1,
          });
          const cylinder = new THREE.Mesh(cylGeo, cylMat);
          cylinder.position.copy(midPoint);

          // Rotate cylinder along bond vector
          const direction = new THREE.Vector3().subVectors(p2, p1).normalize();
          const axis = new THREE.Vector3(0, 1, 0).cross(direction).normalize();
          const angle = Math.acos(new THREE.Vector3(0, 1, 0).dot(direction));
          cylinder.setRotationFromAxisAngle(axis, angle);

          group.add(cylinder);
        }
      });
    }
  }, [conformer, style, showHydrogens]);

  // Mouse Orbit
  const handleMouseDown = (e: React.MouseEvent) => {
    controlsRef.current.isDragging = true;
    controlsRef.current.prevX = e.clientX;
    controlsRef.current.prevY = e.clientY;
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (controlsRef.current.isDragging) {
      const deltaX = e.clientX - controlsRef.current.prevX;
      const deltaY = e.clientY - controlsRef.current.prevY;
      controlsRef.current.rotY -= deltaX * 0.008;
      controlsRef.current.rotX = Math.max(-Math.PI / 2.2, Math.min(Math.PI / 2.2, controlsRef.current.rotX + deltaY * 0.008));
      controlsRef.current.prevX = e.clientX;
      controlsRef.current.prevY = e.clientY;
      return;
    }

    // Atom Raycaster
    if (!mountRef.current || !cameraRef.current || !moleculeGroupRef.current) return;
    const rect = mountRef.current.getBoundingClientRect();
    const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    const mouseY = -((e.clientY - rect.top) / rect.height) * 2 + 1;

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera(new THREE.Vector2(mouseX, mouseY), cameraRef.current);
    const intersects = raycaster.intersectObjects(moleculeGroupRef.current.children);

    if (intersects.length > 0 && intersects[0].object.userData?.atom) {
      setHoveredAtom(intersects[0].object.userData.atom);
      setTooltipPos({ x: e.clientX - rect.left + 15, y: e.clientY - rect.top - 20 });
    } else {
      setHoveredAtom(null);
      setTooltipPos(null);
    }
  };

  const handleMouseUp = () => {
    controlsRef.current.isDragging = false;
  };

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    const zoomDelta = e.deltaY * 0.01;
    controlsRef.current.distance = Math.max(3, Math.min(30, controlsRef.current.distance + zoomDelta));
  };

  const resetCamera = () => {
    controlsRef.current.rotX = 0.2;
    controlsRef.current.rotY = 0.3;
    if (conformer) {
      const maxDim = Math.max(conformer.dimensions.x, conformer.dimensions.y, conformer.dimensions.z, 2.0);
      controlsRef.current.distance = Math.max(7, maxDim * 2.3);
    }
  };

  return (
    <div className="relative w-full h-full flex flex-col bg-slate-950 overflow-hidden select-none border-l border-slate-800">
      {/* 3D Canvas */}
      <div
        ref={mountRef}
        className="w-full flex-1 cursor-grab active:cursor-grabbing"
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onWheel={handleWheel}
      />

      {/* Floating Toolbar */}
      <div className="absolute top-3 left-3 flex flex-wrap gap-1.5 z-10">
        <button
          onClick={resetCamera}
          className="flex items-center gap-1 px-2 py-1 bg-slate-900/90 hover:bg-slate-800 text-slate-300 text-xs font-mono rounded border border-slate-700/70 shadow backdrop-blur transition"
          title="Reset Camera"
        >
          <RotateCcw className="w-3 h-3 text-cyan-400" />
        </button>
        <button
          onClick={() => setAutoRotate(!autoRotate)}
          className={`flex items-center gap-1 px-2 py-1 text-xs font-mono rounded border backdrop-blur transition ${
            autoRotate
              ? 'bg-cyan-500/20 text-cyan-300 border-cyan-500/50'
              : 'bg-slate-900/90 text-slate-400 border-slate-700/70 hover:text-slate-200'
          }`}
          title="Toggle Auto Spin"
        >
          <Compass className="w-3 h-3" />
        </button>
        <button
          onClick={() => setShowHydrogens(!showHydrogens)}
          className={`flex items-center gap-1 px-2 py-1 text-xs font-mono rounded border backdrop-blur transition ${
            showHydrogens
              ? 'bg-indigo-500/20 text-indigo-300 border-indigo-500/50'
              : 'bg-slate-900/90 text-slate-400 border-slate-700/70 hover:text-slate-200'
          }`}
          title="Toggle Hydrogens"
        >
          <span className="font-bold">H</span>
        </button>
      </div>

      {/* Style Selector Tabs */}
      <div className="absolute top-3 right-3 bg-slate-900/90 border border-slate-800 rounded p-1 flex gap-1 z-10 backdrop-blur shadow-md">
        {(['ballStick', 'spacefill', 'sticks', 'wireframe'] as RenderStyle[]).map(s => (
          <button
            key={s}
            onClick={() => setStyle(s)}
            className={`px-2 py-0.5 text-[10px] font-mono rounded transition capitalize ${
              style === s
                ? 'bg-cyan-500 text-slate-950 font-bold'
                : 'text-slate-400 hover:text-slate-200 hover:bg-slate-800'
            }`}
          >
            {s === 'ballStick' ? 'Ball & Stick' : s}
          </button>
        ))}
      </div>

      {/* Conformer Physical Geometry Metrics Footer */}
      {conformer && (
        <div className="px-3 py-2 bg-slate-900/95 border-t border-slate-800 text-xs font-mono text-slate-300 flex items-center justify-between flex-wrap gap-2 z-10">
          <div className="flex items-center gap-3">
            <span className="text-cyan-400 font-bold">{conformer.molecularFormula}</span>
            <span className="text-slate-500">|</span>
            <span className="text-slate-400">
              Atoms: <span className="text-white">{conformer.atoms.length}</span> ({conformer.atoms.filter(a => a.element !== 'H').length} heavy)
            </span>
            <span className="text-slate-500">|</span>
            <span className="text-slate-400">
              Bonds: <span className="text-white">{conformer.bonds.length}</span>
            </span>
          </div>
          <div className="flex items-center gap-3 text-[11px] text-slate-400">
            <span>R_g: <span className="text-emerald-400 font-medium">{conformer.radiusOfGyration} Å</span></span>
            <span>Box: <span className="text-white">{conformer.dimensions.x} × {conformer.dimensions.y} × {conformer.dimensions.z} Å</span></span>
          </div>
        </div>
      )}

      {/* Atom Hover Tooltip */}
      {hoveredAtom && tooltipPos && (
        <div
          className="absolute z-30 pointer-events-none bg-slate-900/95 text-slate-100 p-2 rounded shadow-xl border border-cyan-500/40 text-[11px] font-mono backdrop-blur min-w-[130px]"
          style={{ left: `${tooltipPos.x}px`, top: `${tooltipPos.y}px` }}
        >
          <div className="font-bold text-cyan-300 flex items-center justify-between">
            <span>Atom #{hoveredAtom.id + 1} ({hoveredAtom.symbol})</span>
            {hoveredAtom.aromatic && <span className="text-[9px] text-amber-400">Aromatic</span>}
          </div>
          <div className="text-[10px] text-slate-400 mt-1">
            X: {hoveredAtom.x.toFixed(2)}, Y: {hoveredAtom.y.toFixed(2)}, Z: {hoveredAtom.z.toFixed(2)}
          </div>
        </div>
      )}
    </div>
  );
};
