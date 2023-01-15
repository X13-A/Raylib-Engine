#include "My3DPrimitives.h"
#include "MyConversion.h"
#include <raylib.h>

#pragma region quad

void MyDrawPolygonQuad(Quad quad, Color color)
{
	int numVertex = 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	// BEGINNING OF SPACE TRANSFORMATION INDUCED BY THE LOCAL REFERENCE FRAME
	// methods should be called in this order: rlTranslatef, rlRotatef & rlScalef
	// so that transformations occur in the opposite order: scale, then rotation, then translation
	rlPushMatrix();
	//TRANSLATION
	rlTranslatef(quad.ref.origin.x, quad.ref.origin.y, quad.ref.origin.z);
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(quad.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//SCALING
	rlScalef(quad.extents.x, 1, quad.extents.z);
	// END OF SPACE TRANSFORMATION INDUCED BY THE LOCAL REFERENCE FRAME
	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlVertex3f(1, 0, 1);
	rlVertex3f(1, 0, -1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(1, 0, 1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(-1, 0, 1);
	rlEnd();
	//EVERY rlPushMatrix method call should be followed by a rlPopMatrix method call
	rlPopMatrix();
}

void MyDrawWireframeQuad(Quad quad, Color color)
{
	int numVertex = 10;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	rlPushMatrix();
	rlTranslatef(quad.ref.origin.x, quad.ref.origin.y, quad.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(quad.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(quad.extents.x, 1, quad.extents.z);
	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlVertex3f(1, 0, 1);
	rlVertex3f(1, 0, -1);
	rlVertex3f(1, 0, -1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(1, 0, 1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(-1, 0, 1);
	rlVertex3f(-1, 0, 1);
	rlVertex3f(1, 0, 1);
	rlEnd();
	rlPopMatrix();
}

void MyDrawQuad(Quad quad, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonQuad(quad, polygonColor);
	if (drawWireframe)MyDrawWireframeQuad(quad, wireframeColor);
}

#pragma endregion

#pragma region box

void MyDrawPolygonBox(Box box, Color color)
{
	int numVertex = 36;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(box.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(box.extents.x, box.extents.y, box.extents.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	// Top
	rlVertex3f(-1, 1, -1);
	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, 1, 1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, 1);
	rlVertex3f(1, 1, -1);

	// Bottom
	rlVertex3f(-1, -1, -1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, -1, -1);

	// Front
	rlVertex3f(-1, 1, 1);
	rlVertex3f(-1, -1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, 1, 1);

	// Back
	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, -1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(-1, -1, -1);

	// Left
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, 1);
	rlVertex3f(-1, 1, -1);
	rlVertex3f(-1, -1, -1);
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, -1);

	// Right
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(1, 1, -1);
	rlVertex3f(1, 1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, 1, -1);

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeBox(Box box, Color color)
{
	int numVertex = 36;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(box.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(box.extents.x, box.extents.y, box.extents.z);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	// Front
	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, 1, 1);

	rlVertex3f(1, 1, 1);
	rlVertex3f(1, -1, 1);

	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, -1, 1);

	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, 1);

	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, -1, 1);

	// Top
	rlVertex3f(-1, 1, 1);
	rlVertex3f(-1, 1, -1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, -1);

	rlVertex3f(1, 1, -1);
	rlVertex3f(1, 1, 1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, 1);

	// Back
	rlVertex3f(1, 1, -1);
	rlVertex3f(1, -1, -1);

	rlVertex3f(1, -1, -1);
	rlVertex3f(-1, -1, -1);

	rlVertex3f(-1, -1, -1);
	rlVertex3f(-1, 1, -1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, -1, -1);

	// Bottom
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, -1, -1);

	rlVertex3f(1, -1, 1);
	rlVertex3f(1, -1, -1);

	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, -1, -1);

	// Right
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, 1, -1);

	// Left
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, -1);

	rlEnd();
	rlPopMatrix();
}

void MyDrawBox(Box box, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonBox(box, polygonColor);
	if (drawWireframe) MyDrawWireframeBox(box, wireframeColor);
}

#pragma endregion

#pragma region disk

void MyDrawPolygonDisk(Disk disk, int nSectors, Color color)
{
	int numVertex = nSectors * 3;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(disk.ref.origin.x, disk.ref.origin.y, disk.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(disk.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(disk.radius, 1, disk.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = 0;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(0, 0, 0);
		rlVertex3f(point1.x, point1.y, point1.z);
		rlVertex3f(point2.x, point2.y, point2.z);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeDisk(Disk disk, int nSectors, Color color)
{
	int numVertex = nSectors * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(disk.ref.origin.x, disk.ref.origin.y, disk.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(disk.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(disk.radius, 1, disk.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = 0;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(0, 0, 0);
		rlVertex3f(point1.x, point1.y, point1.z);

		rlVertex3f(point1.x, point1.y, point1.z);
		rlVertex3f(point2.x, point2.y, point2.z);

		rlVertex3f(point2.x, point2.y, point2.z);
		rlVertex3f(0, 0, 0);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawDisk(Disk disk, int nSectors, bool drawPolygon, bool drawWireframe , Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonDisk(disk, nSectors, polygonColor);
	if (drawWireframe) MyDrawWireframeDisk(disk, nSectors, wireframeColor);
}

#pragma endregion

#pragma region sphere

void MyDrawPolygonSphere(Sphere sphere, int nMeridians, int nParallels, Color color)
{
	int numVertex = 6 * nParallels * nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphere.ref.origin.x, sphere.ref.origin.y, sphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphere.radius, sphere.radius, sphere.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / nParallels);
	float meridianAngle = (2 * PI / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
		{
			// Premier triangle
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			// Second triangle
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeSphere(Sphere sphere, int nMeridians, int nParallels, Color color)
{
	int numVertex = 6 * nParallels * nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphere.ref.origin.x, sphere.ref.origin.y, sphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphere.radius, sphere.radius, sphere.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / nParallels);
	float meridianAngle = (2 * PI / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
		{
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawSphere(Sphere sphere, int nMeridians, int nParallels, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonSphere(sphere, nMeridians, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeSphere(sphere, nMeridians, nParallels, wireframeColor);
}

#pragma endregion

#pragma region sphereCorner

void MyDrawPolygonSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, Color color)
{
	int numVertex = 6 * nParallels * nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphereCorner.ref.origin.x, sphereCorner.ref.origin.y, sphereCorner.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphereCorner.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphereCorner.radius, sphereCorner.radius, sphereCorner.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = ((PI / 2) / nParallels);
	float meridianAngle = (2 * (PI / 4) / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
		{
			// Premier triangle
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			//// Second triangle
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, Color color)
{
	int numVertex = (6 * nParallels * nMeridians) + (2 * nMeridians) + (2 * nParallels);
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphereCorner.ref.origin.x, sphereCorner.ref.origin.y, sphereCorner.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphereCorner.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphereCorner.radius, sphereCorner.radius, sphereCorner.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = ((PI / 2) / nParallels);
	float meridianAngle = (2 * (PI / 4) / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
		{
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			if (i == nMeridians - 1) {
				sphPoint.theta = meridianAngle * i;
				sphPoint.phi = parallelAngle * j;
				sphPoint.theta += meridianAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);

				sphPoint.phi += parallelAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);
			}
			if (j == nParallels - 1) {
				sphPoint.theta = meridianAngle * i;
				sphPoint.phi = parallelAngle * j;
				sphPoint.phi += parallelAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);

				sphPoint.theta += meridianAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);
			}
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonSphereCorner(sphereCorner, nMeridians, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeSphereCorner(sphereCorner, nMeridians, nParallels, wireframeColor);
}

#pragma endregion

#pragma region hemisphere

void MyDrawWireframeHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, Color color)
{
	int numVertex = (6 * nParallels * nMeridians) + (2 * nMeridians);
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(hemisphere.ref.origin.x, hemisphere.ref.origin.y, hemisphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(hemisphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(hemisphere.radius, hemisphere.radius, hemisphere.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / 2) / nParallels;
	float meridianAngle = (2 * PI) / nMeridians;

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
		{
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			if (j == nParallels - 1) {
				sphPoint.theta = meridianAngle * i;
				sphPoint.phi = parallelAngle * j;
				sphPoint.phi += parallelAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);

				sphPoint.theta += meridianAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);
			}
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawPolygonHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, Color color)
{
	int numVertex = 6 * nParallels * nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(hemisphere.ref.origin.x, hemisphere.ref.origin.y, hemisphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(hemisphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(hemisphere.radius, hemisphere.radius, hemisphere.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / 2) / nParallels;
	float meridianAngle = (2 * PI) / nMeridians;

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
		{
			// Premier triangle
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			//// Second triangle
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonHemisphere(hemisphere, nMeridians, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeHemisphere(hemisphere, nMeridians, nParallels, wireframeColor);
}

#pragma endregion

#pragma region cylinder

void MyDrawPolygonCylinder(Cylinder cylinder, int nSectors, bool drawCaps, Color color) {
	int numVertex = nSectors * 6;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinder.ref.origin.x, cylinder.ref.origin.y, cylinder.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinder.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinder.radius, cylinder.halfHeight, cylinder.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

		// Premier triangle
		rlVertex3f(point2.x, 1, point2.z);
		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point1.x, -1, point1.z);

		// Second triangle
		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);
		rlVertex3f(point2.x, 1, point2.z);

		// Fermeture du cylindre
		if (drawCaps) {
			rlVertex3f(0, 1, 0);
			rlVertex3f(point1.x, 1, point1.z);
			rlVertex3f(point2.x, 1, point2.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point2.x, -1, point2.z);
			rlVertex3f(point1.x, -1, point1.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeCylinder(Cylinder cylinder, int nSectors, bool drawCaps, Color color) {
	int numVertex = nSectors * 8;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinder.ref.origin.x, cylinder.ref.origin.y, cylinder.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinder.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinder.radius, cylinder.halfHeight, cylinder.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);

		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point1.x, 1, point1.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);

		if (drawCaps) {
			rlVertex3f(0, 1, 0);
			rlVertex3f(point1.x, 1, point1.z);

			rlVertex3f(0, 1, 0);
			rlVertex3f(point2.x, 1, point2.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point1.x, -1, point1.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point2.x, -1, point2.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinder(Cylinder cylinder, int nSectors, bool drawCaps, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonCylinder(cylinder, nSectors, drawCaps, polygonColor);
	if (drawWireframe) MyDrawWireframeCylinder(cylinder, nSectors, drawCaps, wireframeColor);
}

#pragma endregion

#pragma region cylinderCorner

void MyDrawPolygonCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps, Color color) {
	int numVertex = nSectors * 6;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = (PI / 2) / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinderQuarter.ref.origin.x, cylinderQuarter.ref.origin.y, cylinderQuarter.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinderQuarter.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinderQuarter.radius, cylinderQuarter.halfHeight, cylinderQuarter.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	cylPoint.y = -1;

	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

		// Premier triangle
		rlVertex3f(point2.x, 1, point2.z);
		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point1.x, -1, point1.z);

		// Second triangle
		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);
		rlVertex3f(point2.x, 1, point2.z);

		// Fermeture du cylindre
		if (drawCaps) {
			rlVertex3f(0, 1, 0);
			rlVertex3f(point1.x, 1, point1.z);
			rlVertex3f(point2.x, 1, point2.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point2.x, -1, point2.z);
			rlVertex3f(point1.x, -1, point1.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps, Color color) {
	int numVertex = nSectors * 8;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = (PI / 2) / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinderQuarter.ref.origin.x, cylinderQuarter.ref.origin.y, cylinderQuarter.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinderQuarter.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinderQuarter.radius, cylinderQuarter.halfHeight, cylinderQuarter.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	cylPoint.y = -1;
	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);

		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point1.x, 1, point1.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);

		if (drawCaps) {
			rlVertex3f(0, 1, 0);
			rlVertex3f(point1.x, 1, point1.z);

			rlVertex3f(0, 1, 0);
			rlVertex3f(point2.x, 1, point2.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point1.x, -1, point1.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point2.x, -1, point2.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonCylinderQuarter(cylinderQuarter, nSectors, drawCaps, polygonColor);
	if (drawWireframe) MyDrawWireframeCylinderQuarter(cylinderQuarter, nSectors, drawCaps, wireframeColor);
}

#pragma endregion

#pragma region capsule

void MyDrawPolygonCapsule(Capsule capsule, int nSectors, int nParallels, Color color) {
	int numVertexCylinder = nSectors * 6;
	int numVertexHemispheres = 2 * (6 * nParallels * nSectors);
	int numVertex = numVertexCylinder + numVertexHemispheres;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(capsule.ref.origin.x, capsule.ref.origin.y, capsule.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	ReferenceFrame ref = ReferenceFrame(
		{ 0, 0, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Cylinder cylinder = { ref, capsule.halfHeight, capsule.radius };
	MyDrawPolygonCylinder(cylinder, nSectors, false, color);

	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), PI / 2));
	Hemisphere hemisphere = { ref, capsule.radius };
	MyDrawPolygonHemisphere(hemisphere, nSectors, nParallels, color);

	rlRotatef(PI * RAD2DEG, 1, 0, 0);
	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	hemisphere = { ref, capsule.radius };
	MyDrawPolygonHemisphere(hemisphere, nSectors, nParallels, color);


	rlPopMatrix();
}

void MyDrawWireframeCapsule(Capsule capsule, int nSectors, int nParallels, Color color) {
	int numVertexCylinder = nSectors * 6;
	int numVertexHemispheres = 2 * ((6 * nParallels * nSectors) + (2 * nSectors));
	int numVertex = numVertexCylinder + numVertexHemispheres;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();


	rlPushMatrix();
	rlTranslatef(capsule.ref.origin.x, capsule.ref.origin.y, capsule.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	ReferenceFrame ref = ReferenceFrame(
		{ 0, 0, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Cylinder cylinder = { ref, capsule.halfHeight, capsule.radius };
	MyDrawWireframeCylinder(cylinder, nSectors, false, color);

	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), PI / 2));
	Hemisphere hemisphere = { ref, capsule.radius };
	MyDrawWireframeHemisphere(hemisphere, nSectors, nParallels, color);

	rlRotatef(PI * RAD2DEG, 1, 0, 0);
	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	hemisphere = { ref, capsule.radius };
	MyDrawWireframeHemisphere(hemisphere, nSectors, nParallels, color);
	rlRotatef(PI * RAD2DEG, 1, 0, 0);


	rlPopMatrix();
}

void MyDrawCapsule(Capsule capsule, int nSectors, int nParallels, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonCapsule(capsule, nSectors, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeCapsule(capsule, nSectors, nParallels, wireframeColor);
}

#pragma endregion

#pragma region RoundedBox

void MyDrawPolygonRoundedBox(RoundedBox roundedBox, int nSectors, Color color)
{
	int numVertex = 6 * 6; // Quads
	numVertex += 12 * nSectors * 6; // Cylinder quarters
	numVertex += 8 * 6 * nSectors * nSectors; // Sphere corners
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(roundedBox.ref.origin.x, roundedBox.ref.origin.y, roundedBox.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlColor4ub(color.r, color.g, color.b, color.a);

	roundedBox.extents.x /= 2;
	roundedBox.extents.y /= 2;
	roundedBox.extents.z /= 2;
	roundedBox.radius /= 2;

#pragma region Pavé

	ReferenceFrame top_ref = ReferenceFrame(
		{ 0,roundedBox.extents.z + roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Quad top_quad = { top_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(top_quad, color);

	ReferenceFrame bottom_ref = ReferenceFrame(
		{ 0,-roundedBox.extents.z - roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI));
	Quad bottom_quad = { bottom_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(bottom_quad, color);

	ReferenceFrame front_ref = ReferenceFrame(
		{ 0,0,roundedBox.extents.y + roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	Quad front_quad = { front_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(front_quad, color);

	ReferenceFrame back_ref = ReferenceFrame(
		{ 0,0,-roundedBox.extents.y - roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	Quad back_quad = { back_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(back_quad, color);

	ReferenceFrame left_ref = ReferenceFrame(
		{ -roundedBox.extents.x - roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2)));
	Quad left_quad = { left_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(left_quad, color);

	ReferenceFrame right_ref = ReferenceFrame(
		{ roundedBox.extents.x + roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2)));
	Quad right_quad = { right_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(right_quad, color);

#pragma endregion

#pragma region Quart de cylindre vertical

	ReferenceFrame first_vertical_ref(
		{ roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	CylinderQuarter firstVerticalCylinderQuarter = { first_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_vertical_ref(
		{ -roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), 3 * PI / 2));
	CylinderQuarter secondVerticalCylinderQuarter = { second_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_vertical_ref(
		{ -roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	CylinderQuarter thirdVerticalCylinderQuarter = { third_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_vertical_ref(
		{ roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	CylinderQuarter fourthVerticalCylinderQuarter = { fourth_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthVerticalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face supérieure

	ReferenceFrame first_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2));
	CylinderQuarter firstSuperiorHorizontalCylinderQuarter = { first_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondSuperiorHorizontalCylinderQuarter = { second_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_superior_horizontal_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	CylinderQuarter thirdSuperiorHorizontalCylinderQuarter = { third_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_superior_horizontal_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthSuperiorHorizontalCylinderQuarter = { fourth_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthSuperiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face inférieur

	ReferenceFrame first_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2));
	CylinderQuarter firstInferiorHorizontalCylinderQuarter = { first_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondInferiorHorizontalCylinderQuarter = { second_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_inferior_horizontal_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	CylinderQuarter thirdInferiorHorizontalCylinderQuarter = { third_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_inferior_horizontal_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthInferiorHorizontalCylinderQuarter = { fourth_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthInferiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Coin face supérieur

	ReferenceFrame first_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	SphereCorner firstSuperiorCorner = { first_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(firstSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	SphereCorner secondSuperiorCorner = { second_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(secondSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2));
	SphereCorner thirdSuperiorCorner = { third_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(thirdSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	SphereCorner fourthSuperiorCorner = { fourth_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(fourthSuperiorCorner, nSectors, nSectors, color);

#pragma endregion

#pragma region Coin face inférieur

	ReferenceFrame first_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	SphereCorner firstInferiorCorner = { first_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(firstInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner secondInferiorCorner = { second_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(secondInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI)));
	SphereCorner thirdInferiorCorner = { third_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(thirdInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner fourthInferiorCorner = { fourth_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(fourthInferiorCorner, nSectors, nSectors, color);

#pragma endregion

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeRoundedBox(RoundedBox roundedBox, int nSectors, Color color)
{
	int numVertex = 6 * 10; // Quads
	numVertex += 12 * nSectors * 8; // Cylinder quarters
	numVertex += 8 * ((6 * nSectors * nSectors) + (2 * nSectors) + (2 * nSectors)); // Sphere corners

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(roundedBox.ref.origin.x, roundedBox.ref.origin.y, roundedBox.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlColor4ub(color.r, color.g, color.b, color.a);

	roundedBox.extents.x /= 2;
	roundedBox.extents.y /= 2;
	roundedBox.extents.z /= 2;
	roundedBox.radius /= 2;

#pragma region Pavé

	ReferenceFrame top_ref = ReferenceFrame(
		{ 0,roundedBox.extents.z + roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Quad top_quad = { top_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(top_quad, color);

	ReferenceFrame bottom_ref = ReferenceFrame(
		{ 0,-roundedBox.extents.z - roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI));
	Quad bottom_quad = { bottom_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(bottom_quad, color);

	ReferenceFrame front_ref = ReferenceFrame(
		{ 0,0,roundedBox.extents.y + roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	Quad front_quad = { front_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(front_quad, color);

	ReferenceFrame back_ref = ReferenceFrame(
		{ 0,0,-roundedBox.extents.y - roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	Quad back_quad = { back_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(back_quad, color);

	ReferenceFrame left_ref = ReferenceFrame(
		{ -roundedBox.extents.x - roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2)));
	Quad left_quad = { left_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(left_quad, color);

	ReferenceFrame right_ref = ReferenceFrame(
		{ roundedBox.extents.x + roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2)));
	Quad right_quad = { right_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(right_quad, color);

#pragma endregion

#pragma region Quart de cylindre vertical

	ReferenceFrame first_vertical_ref(
		{ roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	CylinderQuarter firstVerticalCylinderQuarter = { first_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_vertical_ref(
		{ -roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), 3 * PI / 2));
	CylinderQuarter secondVerticalCylinderQuarter = { second_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_vertical_ref(
		{ -roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	CylinderQuarter thirdVerticalCylinderQuarter = { third_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_vertical_ref(
		{ roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	CylinderQuarter fourthVerticalCylinderQuarter = { fourth_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthVerticalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face supérieure

	ReferenceFrame first_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2));
	CylinderQuarter firstSuperiorHorizontalCylinderQuarter = { first_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondSuperiorHorizontalCylinderQuarter = { second_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_superior_horizontal_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	CylinderQuarter thirdSuperiorHorizontalCylinderQuarter = { third_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_superior_horizontal_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthSuperiorHorizontalCylinderQuarter = { fourth_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthSuperiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face inférieur

	ReferenceFrame first_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2));
	CylinderQuarter firstInferiorHorizontalCylinderQuarter = { first_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondInferiorHorizontalCylinderQuarter = { second_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_inferior_horizontal_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	CylinderQuarter thirdInferiorHorizontalCylinderQuarter = { third_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_inferior_horizontal_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthInferiorHorizontalCylinderQuarter = { fourth_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthInferiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Coin face supérieur

	ReferenceFrame first_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	SphereCorner firstSuperiorCorner = { first_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(firstSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	SphereCorner secondSuperiorCorner = { second_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(secondSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2));
	SphereCorner thirdSuperiorCorner = { third_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(thirdSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	SphereCorner fourthSuperiorCorner = { fourth_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(fourthSuperiorCorner, nSectors, nSectors, color);

#pragma endregion

#pragma region Coin face inférieur

	ReferenceFrame first_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	SphereCorner firstInferiorCorner = { first_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(firstInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner secondInferiorCorner = { second_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(secondInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI)));
	SphereCorner thirdInferiorCorner = { third_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(thirdInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner fourthInferiorCorner = { fourth_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(fourthInferiorCorner, nSectors, nSectors, color);

#pragma endregion

	rlEnd();
	rlPopMatrix();
}

void MyDrawRoundedBox(RoundedBox roundedBox, int nSectors, bool drawPolygon, bool drawWireframe, Color polygonColor, Color wireframeColor)
{
	if (drawPolygon) MyDrawPolygonRoundedBox(roundedBox, nSectors, polygonColor);
	if (drawWireframe) MyDrawWireframeRoundedBox(roundedBox, nSectors, wireframeColor);
}

#pragma endregion

//#pragma region segment
//
//void MyDrawSegment(Segment segment, Color colorA = RED, Color colorB = GREEN, Color color = BLACK) {
//	DrawLine3D(segment.a, segment.b, color);
//	MyDrawPolygonSphere({ {segment.a,QuaternionIdentity()},.15f }, 16, 8, colorA);
//	MyDrawPolygonSphere({ {segment.b,QuaternionIdentity()},.15f }, 16, 8, colorB);
//
//}
//
//#pragma endregion

