# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@author: P.J.Lawrence
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay


class PlotterClass:
	# data ranges used for ploting
	rangeX = [0, -1]
	rangeY = [0, -1]
	rangeZ = [0, -1]

	def __init__(self):
		self.fig = plt.figure()
		self.plot_area = self.fig.add_subplot(projection='3d')

	def update_ranges(self, lminX, lmaxX, lminY, lmaxY, lminZ, lmaxZ):
		if self.rangeX[1] < self.rangeX[0]:
			self.rangeX = [lminX, lmaxX]
			self.rangeY = [lminY, lmaxY]
			self.rangeZ = [lminZ, lmaxZ]
		else:
			self.rangeX = [min(lminX, self.rangeX[0]), max(lmaxX, self.rangeX[1])]
			self.rangeY = [min(lminY, self.rangeY[0]), max(lmaxY, self.rangeY[1])]
			self.rangeZ = [min(lminZ, self.rangeZ[0]), max(lmaxZ, self.rangeZ[1])]

	def add_element_tri(self, shape):
		"""
		plot the complete object of an entity_shape.geometry
		"""
		faces = shape.faces  # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
		verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
		
		iPoint = 0
		xp = []
		yp = []
		zp = []
		VertexCount = len(verts)
		while iPoint < VertexCount:
			xp.append(verts[iPoint])
			iPoint += 1
			yp.append(verts[iPoint])
			iPoint += 1
			zp.append(verts[iPoint])
			iPoint += 1

		lminX = min(xp)
		lmaxX = max(xp)
		lminY = min(yp)
		lmaxY = max(yp)
		lminZ = min(zp)
		lmaxZ = max(zp)
		
		self.update_ranges(lminX, lmaxX, lminY, lmaxY, lminZ, lmaxZ)
			
		iPoint = 0
		tri_idx = []
		FaceCount = len(faces)
		while iPoint < FaceCount:
			ixp = int(faces[iPoint])
			iPoint += 1
			iyp = int(faces[iPoint])
			iPoint += 1
			izp = int(faces[iPoint])
			iPoint += 1
			tri_idx.append((ixp, iyp, izp))
			
		self.add_faces(xp, yp, zp, tri_idx)

	def add_faces(self, x, y, z, tri_idx):
		self.plot_area.plot_trisurf(x, y, z, triangles=tri_idx)

	def add_trisurf(self, xp, yp, zp, tri_idx):
		self.plot_area.plot_trisurf(xp, yp, zp, triangles=tri_idx)

	def add_tri_delaunay(self, points2D, z):
		mynpPoints = np.array(points2D)
		tri = Delaunay(mynpPoints)
		z_points = [z]*len(points2D)
		self.plot_area.plot_trisurf(mynpPoints[:, 0], mynpPoints[:, 1], z_points, triangles=tri.simplices)

	def show_geom(self):
		distX = self.rangeX[1]-self.rangeX[0]
		distY = self.rangeY[1]-self.rangeY[0]
		distZ = self.rangeZ[1]-self.rangeX[0]
		range_m = (max(distX, distY, distZ))/2
		midX = (self.rangeX[1]+self.rangeX[0])/2
		midY = (self.rangeY[1]+self.rangeY[0])/2
		midZ = (self.rangeZ[1]+self.rangeZ[0])/2

		self.rangeX = [midX - range_m, midX + range_m]
		self.rangeY = [midY - range_m, midY + range_m]
		self.rangeZ = [midZ - range_m, midZ + range_m]
		#print (self.rangeX,self.rangeY,self.rangeZ)
		self.plot_area.set_xlim(self.rangeX)
		self.plot_area.set_ylim(self.rangeY)
		self.plot_area.set_zlim(self.rangeZ)

		plt.show()
