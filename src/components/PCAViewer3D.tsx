/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { MoleculeRecord, ClusterInfo } from '../types/cheminformatics';
import { RotateCcw, ZoomIn, ZoomOut, Compass, Eye, Filter, Sparkles } from 'lucide-react';

interface PCAViewer3DProps {
  molecules: MoleculeRecord[];
  clusters: ClusterInfo[];
  activeMolecule: MoleculeRecord | null;
  onSelectMolecule: (m: MoleculeRecord) => void;
  nearestNeighbors?: Array<{ id: string; name: string; distance: number; approved: boolean }>;
}

export const PCAViewer3D: React.FC<PCAViewer3DProps> = ({
  molecules,
  clusters,
  activeMolecule,
  onSelectMolecule,
  nearestNeighbors = [],
}) => {
  const mountRef = useRef<HTMLDivElement>(null);
  const [hoveredMol, setHoveredMol] = useState<MoleculeRecord | null>(null);
  const [tooltipPos, setTooltipPos] = useState<{ x: number; y: number } | null>(null);
  const [showApprovedOnly, setShowApprovedOnly] = useState(false);
  const [showAxes, setShowAxes] = useState(true);
  const [showNeighborLines, setShowNeighborLines] = useState(true);
  const [autoRotate, setAutoRotate] = useState(false);

  // Three.js internal references
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
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
    rotX: 0.35,
    rotY: 0.75,
    distance: 14,
  });

  const pointsGroupRef = useRef<THREE.Group | null>(null);
  const beaconMeshRef = useRef<THREE.Mesh | null>(null);
  const neighborLinesGroupRef = useRef<THREE.Group | null>(null);

  // Initialize Three.js scene
  useEffect(() => {
    if (!mountRef.current) return;
    const width = mountRef.current.clientWidth;
    const height = mountRef.current.clientHeight;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0b1120); // Deep scientific slate
    sceneRef.current = scene;

    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 100);
    cameraRef.current = camera;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    mountRef.current.innerHTML = '';
    mountRef.current.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Ambient & directional lights
    const ambLight = new THREE.AmbientLight(0xffffff, 0.85);
    scene.add(ambLight);

    const dirLight1 = new THREE.DirectionalLight(0xffffff, 1.2);
    dirLight1.position.set(10, 20, 15);
    scene.add(dirLight1);

    const dirLight2 = new THREE.DirectionalLight(0x60a5fa, 0.6);
    dirLight2.position.set(-15, -10, -10);
    scene.add(dirLight2);

    // Grid Floor
    const gridHelper = new THREE.GridHelper(16, 16, 0x1e293b, 0x0f172a);
    gridHelper.position.y = -3.5;
    scene.add(gridHelper);

    // Groups
    const pointsGroup = new THREE.Group();
    scene.add(pointsGroup);
    pointsGroupRef.current = pointsGroup;

    const neighborLinesGroup = new THREE.Group();
    scene.add(neighborLinesGroup);
    neighborLinesGroupRef.current = neighborLinesGroup;

    // Pulse Beacon for query molecule
    const beaconGeo = new THREE.SphereGeometry(0.35, 24, 24);
    const beaconMat = new THREE.MeshStandardMaterial({
      color: 0x38bdf8,
      emissive: 0x0284c7,
      emissiveIntensity: 0.9,
      roughness: 0.2,
      metalness: 0.8,
    });
    const beaconMesh = new THREE.Mesh(beaconGeo, beaconMat);
    beaconMesh.visible = false;
    scene.add(beaconMesh);
    beaconMeshRef.current = beaconMesh;

    // Halo ring
    const ringGeo = new THREE.RingGeometry(0.45, 0.65, 32);
    const ringMat = new THREE.MeshBasicMaterial({
      color: 0x38bdf8,
      side: THREE.DoubleSide,
      transparent: true,
      opacity: 0.6,
    });
    const ringMesh = new THREE.Mesh(ringGeo, ringMat);
    ringMesh.rotation.x = Math.PI / 2;
    beaconMesh.add(ringMesh);

    // Axes
    const axesGroup = new THREE.Group();
    axesGroup.name = 'axesGroup';
    const axisLen = 4.5;
    // PC1: Red/Cyan
    const dirX = new THREE.Vector3(1, 0, 0);
    axesGroup.add(new THREE.ArrowHelper(dirX, new THREE.Vector3(0, 0, 0), axisLen, 0xef4444, 0.4, 0.2));
    // PC2: Green/Emerald
    const dirY = new THREE.Vector3(0, 1, 0);
    axesGroup.add(new THREE.ArrowHelper(dirY, new THREE.Vector3(0, 0, 0), axisLen, 0x10b981, 0.4, 0.2));
    // PC3: Blue/Violet
    const dirZ = new THREE.Vector3(0, 0, 1);
    axesGroup.add(new THREE.ArrowHelper(dirZ, new THREE.Vector3(0, 0, 0), axisLen, 0x3b82f6, 0.4, 0.2));
    scene.add(axesGroup);

    // Orbit camera animation loop
    let animId: number;
    let clock = new THREE.Clock();

    const animate = () => {
      animId = requestAnimationFrame(animate);
      const elapsedTime = clock.getElapsedTime();

      if (autoRotate) {
        controlsRef.current.rotY += 0.005;
      }

      // Update camera position from spherical coords
      const ctrl = controlsRef.current;
      const x = ctrl.distance * Math.sin(ctrl.rotY) * Math.cos(ctrl.rotX);
      const y = ctrl.distance * Math.sin(ctrl.rotX);
      const z = ctrl.distance * Math.cos(ctrl.rotY) * Math.cos(ctrl.rotX);

      camera.position.set(x, y, z);
      camera.lookAt(0, 0, 0);

      // Animate beacon glow pulse
      if (beaconMeshRef.current && beaconMeshRef.current.visible) {
        const pulse = 1 + 0.25 * Math.sin(elapsedTime * 6);
        beaconMeshRef.current.scale.set(pulse, pulse, pulse);
        ringMesh.rotation.z += 0.02;
      }

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

  // Update Points when molecules, clusters, or filters change
  useEffect(() => {
    if (!pointsGroupRef.current) return;
    const group = pointsGroupRef.current;
    while (group.children.length > 0) {
      const child = group.children[0];
      group.remove(child);
      if (child instanceof THREE.Mesh) {
        child.geometry.dispose();
        if (Array.isArray(child.material)) child.material.forEach(m => m.dispose());
        else child.material.dispose();
      }
    }

    // Centroid colors from cluster info
    const clusterColorMap: Record<number, number> = {};
    clusters.forEach(c => {
      clusterColorMap[c.id] = parseInt(c.color.replace('#', '0x'), 16);
    });

    // Render individual molecule spheres
    molecules.forEach(mol => {
      if (showApprovedOnly && !mol.approved && !mol.isUserQuery) return;

      const [px, py, pz] = mol.pca;
      // Scale radius by SBOF score / drug-score
      const baseRadius = mol.isUserQuery ? 0.35 : 0.12 + Math.min(0.28, (mol.sbofScore || mol.drugScore) * 0.4);
      const sphereGeo = new THREE.SphereGeometry(baseRadius, 16, 16);

      let colorHex = clusterColorMap[mol.cluster] || 0x06b6d4;
      if (mol.isUserQuery) {
        colorHex = 0x38bdf8; // Neon sky blue
      }

      const mat = new THREE.MeshStandardMaterial({
        color: colorHex,
        roughness: 0.35,
        metalness: mol.approved ? 0.7 : 0.2,
        emissive: mol.isUserQuery ? 0x0284c7 : 0x000000,
        emissiveIntensity: mol.isUserQuery ? 0.8 : 0,
      });

      const mesh = new THREE.Mesh(sphereGeo, mat);
      mesh.position.set(px, py, pz);
      mesh.userData = { molecule: mol };
      group.add(mesh);
    });

    // Render Centroid Polyhedra
    clusters.forEach(c => {
      const [cx, cy, cz] = c.centroid;
      const cGeo = new THREE.OctahedronGeometry(0.38, 0);
      const cMat = new THREE.MeshStandardMaterial({
        color: parseInt(c.color.replace('#', '0x'), 16),
        emissive: parseInt(c.color.replace('#', '0x'), 16),
        emissiveIntensity: 0.5,
        wireframe: true,
      });
      const cMesh = new THREE.Mesh(cGeo, cMat);
      cMesh.position.set(cx, cy, cz);
      group.add(cMesh);
    });
  }, [molecules, clusters, showApprovedOnly]);

  // Position Active Molecule Beacon & Neighbor Lines
  useEffect(() => {
    if (!beaconMeshRef.current || !neighborLinesGroupRef.current) return;
    const linesGroup = neighborLinesGroupRef.current;
    while (linesGroup.children.length > 0) {
      const child = linesGroup.children[0];
      linesGroup.remove(child);
      if (child instanceof THREE.Line) {
        child.geometry.dispose();
        (child.material as THREE.Material).dispose();
      }
    }

    if (activeMolecule) {
      const [ax, ay, az] = activeMolecule.pca;
      beaconMeshRef.current.position.set(ax, ay, az);
      beaconMeshRef.current.visible = true;

      // Draw dashed lines to top 3 nearest neighbors
      if (showNeighborLines && nearestNeighbors.length > 0) {
        nearestNeighbors.slice(0, 4).forEach(nn => {
          const targetMol = molecules.find(m => m.id === nn.id);
          if (targetMol) {
            const [tx, ty, tz] = targetMol.pca;
            const points = [new THREE.Vector3(ax, ay, az), new THREE.Vector3(tx, ty, tz)];
            const lineGeo = new THREE.BufferGeometry().setFromPoints(points);
            const lineMat = new THREE.LineDashedMaterial({
              color: 0x38bdf8,
              dashSize: 0.25,
              gapSize: 0.15,
              linewidth: 2,
            });
            const line = new THREE.Line(lineGeo, lineMat);
            line.computeLineDistances();
            linesGroup.add(line);
          }
        });
      }
    } else {
      beaconMeshRef.current.visible = false;
    }
  }, [activeMolecule, molecules, nearestNeighbors, showNeighborLines]);

  // Toggle Axes
  useEffect(() => {
    if (!sceneRef.current) return;
    const axes = sceneRef.current.getObjectByName('axesGroup');
    if (axes) axes.visible = showAxes;
  }, [showAxes]);

  // Mouse & Orbit interaction
  const handleMouseDown = (e: React.MouseEvent) => {
    controlsRef.current.isDragging = true;
    controlsRef.current.prevX = e.clientX;
    controlsRef.current.prevY = e.clientY;
  };

  const handleMouseMove = (e: React.MouseEvent) => {
    if (controlsRef.current.isDragging) {
      const deltaX = e.clientX - controlsRef.current.prevX;
      const deltaY = e.clientY - controlsRef.current.prevY;
      controlsRef.current.rotY -= deltaX * 0.007;
      controlsRef.current.rotX = Math.max(-Math.PI / 2.2, Math.min(Math.PI / 2.2, controlsRef.current.rotX + deltaY * 0.007));
      controlsRef.current.prevX = e.clientX;
      controlsRef.current.prevY = e.clientY;
      return;
    }

    // Raycasting for hover tooltip
    if (!mountRef.current || !cameraRef.current || !pointsGroupRef.current) return;
    const rect = mountRef.current.getBoundingClientRect();
    const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    const mouseY = -((e.clientY - rect.top) / rect.height) * 2 + 1;

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera(new THREE.Vector2(mouseX, mouseY), cameraRef.current);
    const intersects = raycaster.intersectObjects(pointsGroupRef.current.children);

    if (intersects.length > 0 && intersects[0].object.userData?.molecule) {
      setHoveredMol(intersects[0].object.userData.molecule);
      setTooltipPos({ x: e.clientX - rect.left + 15, y: e.clientY - rect.top - 20 });
    } else {
      setHoveredMol(null);
      setTooltipPos(null);
    }
  };

  const handleMouseUp = () => {
    controlsRef.current.isDragging = false;
  };

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    const zoomDelta = e.deltaY * 0.012;
    controlsRef.current.distance = Math.max(4, Math.min(28, controlsRef.current.distance + zoomDelta));
  };

  const handleClick = (e: React.MouseEvent) => {
    if (!mountRef.current || !cameraRef.current || !pointsGroupRef.current) return;
    const rect = mountRef.current.getBoundingClientRect();
    const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    const mouseY = -((e.clientY - rect.top) / rect.height) * 2 + 1;

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera(new THREE.Vector2(mouseX, mouseY), cameraRef.current);
    const intersects = raycaster.intersectObjects(pointsGroupRef.current.children);

    if (intersects.length > 0 && intersects[0].object.userData?.molecule) {
      onSelectMolecule(intersects[0].object.userData.molecule);
    }
  };

  const resetCamera = () => {
    controlsRef.current.rotX = 0.35;
    controlsRef.current.rotY = 0.75;
    controlsRef.current.distance = 14;
  };

  return (
    <div className="relative w-full h-full flex flex-col bg-slate-950 overflow-hidden select-none">
      {/* 3D Canvas Mount */}
      <div
        ref={mountRef}
        className="w-full flex-1 cursor-grab active:cursor-grabbing"
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onWheel={handleWheel}
        onClick={handleClick}
      />

      {/* Floating Viewport Controls */}
      <div className="absolute top-3 left-3 flex flex-wrap gap-2 z-10">
        <button
          onClick={resetCamera}
          className="flex items-center gap-1.5 px-2.5 py-1.5 bg-slate-900/90 hover:bg-slate-800 text-slate-300 text-xs font-mono rounded border border-slate-700/70 shadow-md backdrop-blur transition"
          title="Reset Camera Orientation"
        >
          <RotateCcw className="w-3.5 h-3.5 text-cyan-400" />
          <span>Reset</span>
        </button>
        <button
          onClick={() => setAutoRotate(!autoRotate)}
          className={`flex items-center gap-1.5 px-2.5 py-1.5 text-xs font-mono rounded border backdrop-blur transition ${
            autoRotate
              ? 'bg-cyan-500/20 text-cyan-300 border-cyan-500/50'
              : 'bg-slate-900/90 text-slate-400 border-slate-700/70 hover:text-slate-200'
          }`}
          title="Toggle Auto Rotation"
        >
          <Compass className="w-3.5 h-3.5" />
          <span>Rotate</span>
        </button>
        <button
          onClick={() => setShowApprovedOnly(!showApprovedOnly)}
          className={`flex items-center gap-1.5 px-2.5 py-1.5 text-xs font-mono rounded border backdrop-blur transition ${
            showApprovedOnly
              ? 'bg-emerald-500/20 text-emerald-300 border-emerald-500/50'
              : 'bg-slate-900/90 text-slate-400 border-slate-700/70 hover:text-slate-200'
          }`}
          title="Filter FDA Approved Only"
        >
          <Filter className="w-3.5 h-3.5" />
          <span>FDA Only</span>
        </button>
        <button
          onClick={() => setShowAxes(!showAxes)}
          className={`flex items-center gap-1.5 px-2.5 py-1.5 text-xs font-mono rounded border backdrop-blur transition ${
            showAxes
              ? 'bg-indigo-500/20 text-indigo-300 border-indigo-500/50'
              : 'bg-slate-900/90 text-slate-400 border-slate-700/70 hover:text-slate-200'
          }`}
          title="Toggle PC Axes"
        >
          <Eye className="w-3.5 h-3.5" />
          <span>Axes</span>
        </button>
      </div>

      {/* PCA Loadings Legend */}
      <div className="absolute top-3 right-3 bg-slate-900/90 border border-slate-800 rounded px-3 py-2 text-[11px] font-mono text-slate-300 backdrop-blur z-10 shadow-lg pointer-events-none">
        <div className="font-semibold text-slate-200 mb-1 flex items-center gap-1.5">
          <Sparkles className="w-3 h-3 text-cyan-400" />
          PCA Space (78.4% Var)
        </div>
        <div className="flex items-center gap-2">
          <span className="w-2.5 h-2.5 rounded-full bg-red-500 inline-block" />
          <span className="text-slate-400">PC1 (41.2%):</span> MW, RTB, PSA
        </div>
        <div className="flex items-center gap-2">
          <span className="w-2.5 h-2.5 rounded-full bg-emerald-500 inline-block" />
          <span className="text-slate-400">PC2 (23.8%):</span> logP vs Polarity
        </div>
        <div className="flex items-center gap-2">
          <span className="w-2.5 h-2.5 rounded-full bg-blue-500 inline-block" />
          <span className="text-slate-400">PC3 (13.4%):</span> QED, HBD Donors
        </div>
      </div>

      {/* Cluster Legend Footer */}
      <div className="px-4 py-2 bg-slate-900/90 border-t border-slate-800/80 flex items-center justify-between text-xs font-mono text-slate-300 z-10">
        <div className="flex items-center gap-4 flex-wrap">
          <span className="text-slate-400 text-[11px] uppercase tracking-wider font-semibold">Clusters:</span>
          {clusters.map(c => (
            <div key={c.id} className="flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ backgroundColor: c.color }} />
              <span className="font-medium text-slate-200">{c.name}</span>
              <span className="text-slate-400 text-[10px]">({c.count} cpds, ~{c.meanMw} Da)</span>
            </div>
          ))}
        </div>
        <div className="text-[11px] text-slate-400 hidden sm:block">
          Point Size ∝ SBOF Learned Probability (β=4.0)
        </div>
      </div>

      {/* Hover Tooltip */}
      {hoveredMol && tooltipPos && (
        <div
          className="absolute z-30 pointer-events-none bg-slate-900/95 text-slate-100 p-2.5 rounded shadow-xl border border-cyan-500/40 text-xs font-mono min-w-[210px] backdrop-blur"
          style={{ left: `${tooltipPos.x}px`, top: `${tooltipPos.y}px` }}
        >
          <div className="flex items-center justify-between font-bold text-cyan-300 border-b border-slate-800 pb-1 mb-1.5">
            <span>{hoveredMol.name}</span>
            {hoveredMol.approved && (
              <span className="text-[9px] px-1 py-0.5 bg-emerald-500/20 text-emerald-300 rounded border border-emerald-500/40">
                FDA
              </span>
            )}
          </div>
          <div className="text-[10px] text-slate-400 mb-1">{hoveredMol.id}</div>
          <div className="grid grid-cols-2 gap-x-2 gap-y-0.5 text-[11px] text-slate-300">
            <div>MW: <span className="text-white">{hoveredMol.descriptors.mw}</span></div>
            <div>logP: <span className="text-white">{hoveredMol.descriptors.logP}</span></div>
            <div>QED: <span className="text-cyan-400 font-semibold">{hoveredMol.descriptors.qed}</span></div>
            <div>SBOF: <span className="text-emerald-400 font-semibold">{((hoveredMol.sbofScore || hoveredMol.drugScore) * 100).toFixed(2)}%</span></div>
          </div>
          <div className="mt-1.5 pt-1 border-t border-slate-800 text-[10px] text-slate-400 flex justify-between">
            <span>Cluster: #{hoveredMol.cluster + 1}</span>
            <span className="text-cyan-400">Click to Inspect</span>
          </div>
        </div>
      )}
    </div>
  );
};
