#pragma once

#include <vector>

struct Triangulation
{
	struct Node
	{
		double ro;
		double z;
		int global_ID;
		bool teta_boundery = false;
		bool q_boundery = false;
	};

	struct Element
	{
		int element_ID;
		int nodes_global_ID[3];
		static int nodes_local_ID[3];
		Node nodes_coordinates[3];
		int number_q_boundery_nodes = 0;

		int getGlobalID(int localID) const { return nodes_global_ID[localID]; };
	};

	static Triangulation CreateTriangulation(double ro_begin, double ro_final,
											 double z_begin, double z_final,
		                              int node_number_ro, int node_number_z);

	int GetGlobalID(const Element& element, int localID) const;

	std::vector<Node> nodes;
	std::vector<Element> elements;
	double ro_begin;
	double ro_final;
	double z_begin;
	double z_final;
};