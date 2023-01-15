#pragma once

#include "rlgl.h"
#include "MyStruct.h"

void MyDrawPolygonQuad(Quad quad, Color color);
void MyDrawWireframeQuad(Quad quad, Color color);
void MyDrawQuad(Quad quad, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonDisk(Disk disk, int nSectors, Color color);
void MyDrawWireframeDisk(Disk disk, int nSectors, Color color);
void MyDrawDisk(Disk disk, int nSectors, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonBox(Box box, Color color);
void MyDrawWireframeBox(Box box, Color color);
void MyDrawBox(Box box, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonSphere(Sphere sphere, int nMeridians, int nParallels, Color color);
void MyDrawWireframeSphere(Sphere sphere, int nMeridians, int nParallels, Color color);
void MyDrawSphere(Sphere sphere, int nMeridians, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, Color color = LIGHTGRAY);
void MyDrawWireframeSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, Color color = DARKGRAY);
void MyDrawSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonCylinder(Cylinder cylinder, int nSectors, bool drawCaps, Color color);
void MyDrawWireframeCylinder(Cylinder cylinder, int nSectors, bool drawCaps, Color color);
void MyDrawCylinder(Cylinder cylinder, int nSectors, bool drawCaps = false, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonCapsule(Capsule capsule, int nSectors, int nParallels, Color color);
void MyDrawWireframeCapsule(Capsule capsule, int nSectors, int nParallels, Color color);
void MyDrawCapsule(Capsule capsule, int nSectors, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonRoundedBox(RoundedBox roundedBox, int nSectors, Color color);
void MyDrawWireframeRoundedBox(RoundedBox roundedBox, int nSectors, Color color);
void MyDrawRoundedBox(RoundedBox roundedBox, int nSectors, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawWireframeHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, Color color = DARKGRAY);
void MyDrawPolygonHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, Color color = DARKGRAY);
void MyDrawHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);

void MyDrawPolygonCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps = false, Color color = LIGHTGRAY);
void MyDrawWireframeCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps = false, Color color = DARKGRAY);
void MyDrawCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps = false, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY);
