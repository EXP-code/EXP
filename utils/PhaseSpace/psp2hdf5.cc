#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include "PSP.H"
#include <H5Cpp.h>
#include "cxxopts.H"		// Option parsing
#include "libvars.H"		// EXP library globals

/**
   Create a Gadget2 HDF5 file from PSP.  This should be easy to
   generalize to other HDF5 phase-space variants, most of which are
   Gadget2 like.
 */
int
main(int argc, char **argv)
{
  std::string pspfile, hdf5file, dir;
  std::vector<std::string> partNames(6);
  std::vector<double>      massValue(6, 0.0);
  std::vector<int>         numbValue(6, 0);
  std::vector<bool>        massFixed(6, true);
  bool verbose = false;


  cxxopts::Options options(argv[0], "\nMake a Gadget2-style HDF5 file from a PSP file and a Gadget template file.\nNo cosmological parameters are set.  No subgrid parameters will be set.\nUse the numerical flags to assign component names to Gadget particle types.\n");

  options.add_options()
    ("h,help", "print this help message")
    ("v,verbose", "verbose debugging output")
    ("haloMM", "assign multimass halo array")
    ("bulgeMM", "assign multimass bulge array")
    ("f,infile", "the PSP file",
     cxxopts::value<std::string>(pspfile))
    ("d,dir", "the PSP data directory",
     cxxopts::value<std::string>(dir))
    ("o,outfile", "the hdf5 file",
     cxxopts::value<std::string>(hdf5file)->default_value("new.hdf5"))
    ("0,gas",    "PSP component name for gas",
     cxxopts::value<std::string>(partNames[0]))
    ("1,halo",   "PSP component name for halo",
     cxxopts::value<std::string>(partNames[1]))
    ("2,disk",   "PSP component name for disk",
     cxxopts::value<std::string>(partNames[2]))
    ("3,bulge",  "PSP component name for bulge",
     cxxopts::value<std::string>(partNames[3]))
    ("4,stars",  "PSP component name for stars",
     cxxopts::value<std::string>(partNames[4]))
    ("5,bndry",  "PSP component name for boundary particles (not used)",
     cxxopts::value<std::string>(partNames[5]))
    ;

  auto vm = options.parse(argc, argv);

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  if (vm.count("verbose")) verbose = true;
  
  if (vm.count("infile") == 0) {
    std::cout << "You need to specify the input PSP file" << std::endl
	      << options.help() << std::endl;
    exit(-1);
  }

  std::map<int, std::string> partMap;
  if (vm.count("gas"  )) partMap[0] = partNames[0];
  if (vm.count("halo" )) partMap[1] = partNames[1];
  if (vm.count("disk" )) partMap[2] = partNames[2];
  if (vm.count("bulge")) partMap[3] = partNames[3];
  if (vm.count("stars")) partMap[4] = partNames[4];
  if (vm.count("bndry")) partMap[5] = partNames[5];

  if (vm.count("haloMM"))  massFixed[1] = false;
  if (vm.count("bulgeMM")) massFixed[3] = false;

  PSPptr psp = PSP::getPSP(pspfile, dir, verbose);

  // Now write a summary
  // -------------------
  if (verbose) {
    psp->PrintSummary(cerr);
    cerr << "\nPSP file <" << pspfile << "> has time <" 
	 << psp->CurrentTime() << ">\n";
  }

				// Dump ascii for each component
				// -----------------------------
  PSPstanza *stanza;
  SParticle* part;

  for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {

    for (auto k : partMap) {
      if (k.second.find(stanza->name) == 0) {
	numbValue[k.first] = stanza->comp.nbod;
	if (massFixed[k.first]) massValue[k.first] = psp->GetParticle()->mass();
	if (verbose) {
	  std::cout << "Found stanza name for Type " << k.first
		    << " with " << numbValue[k.first];
	  if (massFixed[k.first])
	    std::cout << " particles with mass=" << massValue[k.first]
		      << std::endl;
	  else
	    std::cout << " particles with variable mass" << std::endl;
	}
      }
    }
  }
  
  // Try block to detect exceptions raised by any of the calls inside it
  //
  try {
    // Turn off the auto-printing when failure occurs so that we can
    // handle the exceptions
    //
    H5::Exception::dontPrint();
    
    // Create the named file, truncating the existing one if any,
    // using default create and access property lists.
    //
    auto file = std::make_shared<H5::H5File>(hdf5file, H5F_ACC_TRUNC);
    
    // Create the Header group
    //
    auto group = std::make_shared<H5::Group>(file->createGroup("/Header"));
    
    // Create the BoxSize attribute
    {
      double value = 0.0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("BoxSize",
						 H5::PredType::NATIVE_DOUBLE,
						 dspace);
      att.write(H5::PredType::NATIVE_DOUBLE, &value);
    }

    // Create the Flag_Cooling attribute
    {
      int value = 0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Flag_Cooling",
						 H5::PredType::NATIVE_INT,
						 dspace);
      att.write(H5::PredType::NATIVE_INT, &value);
    }

    // Create the Flag_Entropy_ICs attribute
    {
      hsize_t dim[] = {6};
      std::vector<unsigned> data(6, 0);
      
      // The length of dim (which is one here)
      int rank = sizeof(dim) / sizeof(hsize_t);

      // preparation of a dataset and a file.
      H5::DataSpace space(rank, dim);
      H5::Attribute att = group->createAttribute("Flag_Entropy_ICs",
						 H5::PredType::NATIVE_UINT,
						 space);
      att.write(H5::PredType::NATIVE_UINT, &data[0]);
    }

    // Create the feedback flag
    {
      int value = 0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Flag_Feedback",
						 H5::PredType::NATIVE_INT,
						 dspace);
      att.write(H5::PredType::NATIVE_INT, &value);
    }

    // Create the metals flag
    {
      int value = 0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Flag_Metals",
						 H5::PredType::NATIVE_INT,
						 dspace);
      att.write(H5::PredType::NATIVE_INT, &value);
    }

    // Create the star-formation flag
    {
      int value = 0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Flag_Sfr",
						 H5::PredType::NATIVE_INT,
						 dspace);
      att.write(H5::PredType::NATIVE_INT, &value);
    }

    // Create the stellar-age flag
    {
      int value = 0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Flag_StellarAge",
						 H5::PredType::NATIVE_INT,
						 dspace);
      att.write(H5::PredType::NATIVE_INT, &value);
    }

    // Create the Hubble parameter
    {
      double value = 0.0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("HubbleParam",
					     H5::PredType::NATIVE_DOUBLE,
					     dspace);
      att.write(H5::PredType::NATIVE_DOUBLE, &value);
    }

    // Create the mass table
    {
      hsize_t dim[] = {6};
      
      // The length of dim (which is one here)
      int rank = sizeof(dim) / sizeof(hsize_t);

      // preparation of a dataset and a file.
      H5::DataSpace space(rank, dim);
      H5::Attribute att = group->createAttribute("MassTable",
						 H5::PredType::NATIVE_DOUBLE,
						 space);
      att.write(H5::PredType::NATIVE_DOUBLE, &massValue[0]);
    }


    // Create the stellar-age flag
    {
      int value = 1;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("NumFilesPerSnapShot",
						 H5::PredType::NATIVE_INT,
						 dspace);
      att.write(H5::PredType::NATIVE_INT, &value);
    }


    // Create the particle number array
    {
      hsize_t dim[] = {6};
      
      // The length of dim (which is one here)
      int rank = sizeof(dim) / sizeof(hsize_t);

      // preparation of a dataset and a file.
      H5::DataSpace space(rank, dim);
      H5::Attribute att = group->createAttribute("NumPart_ThisFile",
						 H5::PredType::NATIVE_INT,
						 space);
      att.write(H5::PredType::NATIVE_INT, &numbValue[0]);
    }

    // Create the total particle number array
    {
      hsize_t dim[] = {6};
      std::vector<unsigned> data;
      for (auto v : numbValue) data.push_back(v);
      
      // The length of dim (which is one here)
      int rank = sizeof(dim) / sizeof(hsize_t);

      // preparation of a dataset and a file.
      H5::DataSpace space(rank, dim);
      H5::Attribute att = group->createAttribute("NumPart_Total",
						 H5::PredType::NATIVE_UINT,
						 space);
      att.write(H5::PredType::NATIVE_UINT, &data[0]);
    }

    // Create the high-word array
    {
      hsize_t dim[] = {6};
      std::vector<unsigned> data(6, 0);
      
      // The length of dim (which is one here)
      int rank = sizeof(dim) / sizeof(hsize_t);

      // preparation of a dataset and a file.
      H5::DataSpace space(rank, dim);
      H5::Attribute att = group->createAttribute("NumPart_Total_HighWord",
						 H5::PredType::NATIVE_UINT,
						 space);
      att.write(H5::PredType::NATIVE_UINT, &data[0]);
    }


    // Create the Omega value
    {
      double value = 0.0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Omega0",
						 H5::PredType::NATIVE_DOUBLE,
						 dspace);
      att.write(H5::PredType::NATIVE_DOUBLE, &value);
    }

    // Create the OmegaLambda value
    {
      double value = 0.0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("OmegaLambda",
						 H5::PredType::NATIVE_DOUBLE,
						 dspace);
      att.write(H5::PredType::NATIVE_DOUBLE, &value);
    }


    // Create the Redshift value
    {
      double value = 0.0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Redshift",
						 H5::PredType::NATIVE_DOUBLE,
						 dspace);
      att.write(H5::PredType::NATIVE_DOUBLE, &value);
    }

    // Create the time value
    {
      double value = 0.0;
      H5::DataSpace dspace(H5S_SCALAR);
      H5::Attribute att = group->createAttribute("Time",
						 H5::PredType::NATIVE_DOUBLE,
						 dspace);
      att.write(H5::PredType::NATIVE_DOUBLE, &value);
    }

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    

      for (auto k : partMap) {

	if (k.second.find(stanza->name) == 0) {

	  // Create coordinates and indices
	  //
	  std::vector<float>    posv(3*stanza->comp.nbod);
	  std::vector<float>    velv(3*stanza->comp.nbod);
	  std::vector<unsigned> indx(  stanza->comp.nbod);
	  std::vector<float>    masv;

	  if (not massFixed[k.first]) masv.resize(stanza->comp.nbod);

	  unsigned cntr=0;
	  for (auto part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	    if (stanza->index_size) indx[cntr] = part->indx();
	    else indx[cntr] = cntr + 1;

	    for (int i=0; i<3; i++)
	      posv[3*cntr + i] = part->pos(i);

	    for (int i=0; i<3; i++)
	      velv[3*cntr + i] = part->vel(i);
	    
	    if (not massFixed[k.first]) masv[cntr] = part->mass();

	    cntr++;
	  }
      
	  // Create the PartType group
	  //
	  std::ostringstream str;
	  str << "/PartType" << k.first;
	  auto partgrp = std::make_shared<H5::Group>(file->createGroup(str.str()));
    
	  hsize_t dim2[] = {static_cast<hsize_t>(stanza->comp.nbod), 3};
      
	  // The length of dim (which is one here)
	  int rank = sizeof(dim2) / sizeof(hsize_t);
	  
	  // preparation of a dataset and a file.
	  H5::DataSpace space(rank, dim2);
	  
	  auto dataset1 = std::make_shared<H5::DataSet>
	    (partgrp->createDataSet("Coordinates",
				    H5::PredType::NATIVE_FLOAT,
				    space));
	  
	  dataset1->write(&posv[0], H5::PredType::NATIVE_FLOAT);

	  auto dataset2 = std::make_shared<H5::DataSet>
	    (partgrp->createDataSet("Velocities",
				    H5::PredType::NATIVE_FLOAT,
				    space));
	  
	  dataset2->write(&velv[0], H5::PredType::NATIVE_FLOAT);


	  hsize_t dim1[] = {static_cast<hsize_t>(stanza->comp.nbod)};
      
	  // The length of dim (which is one here)
	  rank = sizeof(dim1) / sizeof(hsize_t);
	  
	  // preparation of a dataset and a file.
	  H5::DataSpace space1(rank, dim1);

	  auto dataset3 = std::make_shared<H5::DataSet>
	    (partgrp->createDataSet("ParticleIDs",
				    H5::PredType::NATIVE_UINT,
				    space1));
	  
	  dataset3->write(&indx[0], H5::PredType::NATIVE_UINT);

	  if (masv.size()) {
	    auto dataset4 = std::make_shared<H5::DataSet>
	    (partgrp->createDataSet("Masses",
				    H5::PredType::NATIVE_FLOAT,
				    space1));
	    dataset4->write(&masv[0], H5::PredType::NATIVE_FLOAT);
	  }

	  break;
	}
	// End of component-found block
      }
      // End of component-search loop
    }
    // End stanza component loop
  }
  // End of try block

  // catch failure caused by the H5File operations
  catch (H5::FileIException error) {
    error.printErrorStack();
    return -1;
  }
  
  // catch failure caused by the DataSet operations
  catch (H5::DataSetIException error) {
    error.printErrorStack();
    return -1;
  }
  
  // catch failure caused by the DataSpace operations
  catch (H5::DataSpaceIException error) {
    error.printErrorStack();
    return -1;
  }
  
  // catch failure caused by the Attribute operations
  catch (H5::AttributeIException error) {
    error.printErrorStack();
    return -1;
  }

  // DONE!
  //
  return 0;
}
