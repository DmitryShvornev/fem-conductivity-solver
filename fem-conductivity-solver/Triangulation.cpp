#include "Triangulation.h"

Triangulation Triangulation::CreateTriangulation(
	double roBegin, double roFinal,
	double zBegin, double zFinal,
	int nodeNumberRo, int nodeNumberZ)
{
	Triangulation triangulation;
	triangulation.ro_begin = roBegin;
	triangulation.ro_final = roFinal;
	triangulation.z_begin = zBegin;
	triangulation.z_final = zFinal;

	double ro_step = (roFinal - roBegin) / (nodeNumberRo - 1.);
	double z_step = (zFinal - zBegin) / (nodeNumberZ - 1.);

	int ID_counter = 0; 
	for (int ro_step_number = 0; ro_step_number < nodeNumberRo; ro_step_number++)
	{
		for (int z_step_number = 0; z_step_number < nodeNumberZ; z_step_number++)
		{
			Node new_node;
			new_node.global_ID = ID_counter;
			ID_counter++;
			new_node.ro = roBegin + ro_step*ro_step_number;
			new_node.z = zBegin + z_step*z_step_number;
			if ((ro_step_number == 0) || (ro_step_number == (nodeNumberRo - 1)))
				new_node.q_boundery = true;
			if ((z_step_number == 0) || (z_step_number == (nodeNumberZ - 1)))
				new_node.teta_boundery = true;
			triangulation.nodes.push_back(new_node);
		}
	}

	int element_counter = 0;
	for (int ro_node_number = 1; ro_node_number < nodeNumberRo; ro_node_number++)
	{
		for (int z_node_number = 0; z_node_number < (nodeNumberZ - 1); z_node_number++)
		{
			Element new_element1;

			new_element1.element_ID = element_counter;
			element_counter++;
			new_element1.nodes_global_ID[0] = triangulation.nodes.at(
				ro_node_number * nodeNumberZ + z_node_number).global_ID;
			new_element1.nodes_global_ID[1] = triangulation.nodes.at(
				ro_node_number * nodeNumberZ + (z_node_number + 1)).global_ID;
			new_element1.nodes_global_ID[2] = triangulation.nodes.at(
				(ro_node_number - 1) * nodeNumberZ + (z_node_number + 1)).global_ID;
			new_element1.nodes_coordinates[0] 
				= triangulation.nodes.at(new_element1.nodes_global_ID[0]);
			new_element1.nodes_coordinates[1]
				= triangulation.nodes.at(new_element1.nodes_global_ID[1]);
			new_element1.nodes_coordinates[2] 
				= triangulation.nodes.at(new_element1.nodes_global_ID[2]);
			if (triangulation.nodes.at(
				ro_node_number * nodeNumberZ + z_node_number).q_boundery == true)
				new_element1.number_q_boundery_nodes++;
			if (triangulation.nodes.at(
				ro_node_number * nodeNumberZ + (z_node_number + 1)).q_boundery == true)
				new_element1.number_q_boundery_nodes++;
			if (triangulation.nodes.at(
				(ro_node_number - 1) * nodeNumberZ + (z_node_number + 1)).q_boundery == true)
				new_element1.number_q_boundery_nodes++;
			triangulation.elements.push_back(new_element1);

			Element new_element2;
			new_element2.element_ID = element_counter;
			element_counter++;
			new_element2.nodes_global_ID[0] = triangulation.nodes.at(
				ro_node_number * nodeNumberZ + z_node_number).global_ID;
			new_element2.nodes_global_ID[1] = triangulation.nodes.at(
				(ro_node_number - 1) * nodeNumberZ + (z_node_number + 1)).global_ID;
			new_element2.nodes_global_ID[2] = triangulation.nodes.at(
				(ro_node_number - 1) * nodeNumberZ + z_node_number).global_ID;
			new_element2.nodes_coordinates[0] 
				= triangulation.nodes.at(new_element2.nodes_global_ID[0]);
			new_element2.nodes_coordinates[1] 
				= triangulation.nodes.at(new_element2.nodes_global_ID[1]);
			new_element2.nodes_coordinates[2] 
				= triangulation.nodes.at(new_element2.nodes_global_ID[2]);
			if (triangulation.nodes.at(
				ro_node_number * nodeNumberZ + z_node_number).q_boundery == true)
				new_element2.number_q_boundery_nodes++;
			if (triangulation.nodes.at(
				(ro_node_number - 1) * nodeNumberZ + (z_node_number + 1)).q_boundery == true)
				new_element2.number_q_boundery_nodes++;
			if (triangulation.nodes.at(
				(ro_node_number - 1) * nodeNumberZ + z_node_number).q_boundery == true)
				new_element2.number_q_boundery_nodes++;
			triangulation.elements.push_back(new_element2);
		}
	}
	return triangulation;
}

int Triangulation::GetGlobalID(const Element& element, int localID) const
{
	return element.nodes_global_ID[localID];
}

int Triangulation::Element::nodes_local_ID[] = { 0, 1, 2 };