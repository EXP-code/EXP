#include <writePVD.H>
#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

void writePVD(const std::string& filename,
	      const std::vector<double>& times,
	      const std::vector<std::string>& files)
{
  // Sanity check
  //
  if (times.size() != files.size()) {
    std::cerr << "Mismatch in file and time arrays" << std::endl;
    exit(-3);
  }

  // Make the top-level XML file
  //
  rapidxml::xml_document<> doc;
  rapidxml::xml_node<>* decl = doc.allocate_node(rapidxml::node_declaration);
  decl->append_attribute(doc.allocate_attribute("version", "1.0"));
  decl->append_attribute(doc.allocate_attribute("encoding", "UTF-8"));
  doc.append_node(decl);

  // Make the main VTKFile node
  //
  rapidxml::xml_node<>* root = doc.allocate_node(rapidxml::node_element, "VTKFile");
  root->append_attribute(doc.allocate_attribute("type", "Collection"));
  root->append_attribute(doc.allocate_attribute("version", "0.1"));
  root->append_attribute(doc.allocate_attribute("byte_order", "LittleEndian"));
  root->append_attribute(doc.allocate_attribute("compressor", "vtkZLibDataCompressor"));
  doc.append_node(root);

  // Make the Collection
  //
  rapidxml::xml_node<>* coll = doc.allocate_node(rapidxml::node_element, "Collection");

  root->append_node(coll);	// Add the Collection to the root node

  for (size_t i=0; i<times.size(); i++) {
    rapidxml::xml_node<>* child = doc.allocate_node(rapidxml::node_element, "DataSet");

    child->append_attribute(doc.allocate_attribute("timestep", std::to_string(times[i]).c_str()));
    child->append_attribute(doc.allocate_attribute("part", "0"));
    child->append_attribute(doc.allocate_attribute("file", files[i].c_str()));
  }

  // Write the property tree to the XML file.
  //
  std::ofstream out(filename);
  if (out) {
    out << doc;
    std::cout << "Wrote PVD file <" << filename.c_str() << "> "
	      << " with " << times.size() << " data sets." << std::endl;
  } else {
    std::cerr << "writePVD: error opening xml file <" << filename << ">" << std::endl;
  }
}
