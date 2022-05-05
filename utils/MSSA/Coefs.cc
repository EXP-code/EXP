#include "Coefs.H"
#include <iostream>
#include <cstring>
#include <cmath>

bool CylCoefs::read(std::istream& in, bool verbose)
{
  if (in.eof()) return false;

  // iostream exception handling
  //
  in.exceptions ( std::istream::failbit | std::istream::badbit );

  // Save initial stream position
  //
  auto curpos = in.tellg();

  // Attempt to read coefficient magic number
  //
  const unsigned int cmagic = 0xc0a57a3;
  unsigned int tmagic;

  try {

    in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

    if (tmagic == cmagic) {

      // YAML size
      //
      unsigned ssize;
      in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

      // Make and read char buffer
      //
      auto buf = std::make_unique<char[]>(ssize+1);
      in.read(buf.get(), ssize);
      buf[ssize] = 0;		// Null terminate

      YAML::Node node = YAML::Load(buf.get());

      // Get parameters
      //
      time = node["time"].as<double>();
      nmax = node["nmax"].as<int>();
      mmax = node["mmax"].as<int>();

      if (verbose)
	std::cerr << "New header: T=" << time << " nmax=" << nmax
		  << " mmax=" << mmax << std::endl;
    } else {

      // Rewind file
      //
      in.clear();
      in.seekg(curpos);

      CylCoefHeader header;
      in.read((char *)&header, sizeof(CylCoefHeader));

      time = header.time;
      nmax = header.nmax;
      mmax = header.mmax;

      if (verbose)
	std::cerr << "Old header: T=" << time << " nmax=" << nmax
		  << " mmax=" << mmax << std::endl;
    }

    cos_c.resize(mmax+1);
    sin_c.resize(mmax+1);

    for (int mm=0; mm<=mmax; mm++) {

      cos_c[mm].resize(nmax);
      in.read((char *)&cos_c[mm][0], sizeof(double)*nmax);

      if (mm) {
	sin_c[mm].resize(nmax);
	in.read((char *)&sin_c[mm][0], sizeof(double)*nmax);
      }
    }
  }
  catch (std::istream::failure e) {
    if (not in.eof())
      std::cerr << "Coefficient read FAILED at T=" << time
		<< ": " << e.what() << std::endl;
    return false;
  }

  if (verbose)
    std::cerr << "Coefficients successfully read at T=" << time << std::endl;

  return true;
}


void CylCoefs::write(std::ostream& out)
{
  // Coefficient magic number
  const unsigned int cmagic = 0xc0a57a3;

  // Will contain the YAML header
  YAML::Node node;

  node["time"  ] = time;
  node["mmax"  ] = mmax;
  node["nmax"  ] = nmax;

  // Serialize the node
  //
  YAML::Emitter y; y << node;

  // Get the size of the string
  //
  unsigned int hsize = strlen(y.c_str());

  // Write magic #
  //
  out.write(reinterpret_cast<const char *>(&cmagic),   sizeof(unsigned int));

  // Write YAML string size
  //
  out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));

  // Write YAML string
  //
  out.write(reinterpret_cast<const char *>(y.c_str()), hsize);

  // Write the coefficients
  //
  for (int mm=0; mm<=mmax; mm++) {

    for (int j=0; j<nmax; j++)
      out.write((const char *)&cos_c[mm][j], sizeof(double));

    if (mm) {

      for (int j=0; j<nmax; j++)
	out.write((const char *)&sin_c[mm][j], sizeof(double));
    }
  }
}


bool SphCoefs::read(std::istream& in)
{
  if (in.eof()) return false;

  // iostream exception handling
  //
  in.exceptions ( std::istream::failbit | std::istream::badbit );

  // Save initial stream position
  //
  auto curpos = in.tellg();

  // Coefficient magic number
  //
  const unsigned int cmagic = 0xc0a57a2;

  // Try to read magic #
  //
  unsigned int tmagic;

  try {

    in.read(reinterpret_cast<char *>(&tmagic), sizeof(unsigned int));

    // Found new-style coefficient file
    //
    if (cmagic == tmagic) {

      // Read YAML string size
      //
      unsigned int hsize;
      in.read(reinterpret_cast<char *>(&hsize), sizeof(unsigned int));

      // Create buffer
      //
      auto buf = std::make_unique<char[]>(hsize+1);

      // Read YAML string
      //
      in.read(buf.get(), hsize);
      buf[hsize] = 0;		// Null terminate

      YAML::Node node = YAML::Load(buf.get());

      // Get parameters
      //
      header.Lmax  = node["lmax"  ].as<int>();
      header.nmax  = node["nmax"  ].as<int>();
      header.tnow  = node["time"  ].as<double>();
      header.scale = node["scale" ].as<double>();
      normed       = node["normed"].as<bool>();

      std::fill(header.id, header.id+64, 0);
      std::string ID = node["id"].as<std::string>();
      strncpy(header.id, ID.c_str(), std::min<int>(64, ID.size()));

      coefs.resize((header.Lmax+1)*(header.Lmax+1));
      for (auto & v : coefs) v.resize(header.nmax);

      // Read the coefficients
      //
      for (int ir=0; ir<header.nmax; ir++) {
	for (int l=0, loffset=0; l<=header.Lmax; loffset+=(2*l+1), l++) {
	  double fac1 = (2.0*l+1.0)/(4.0*M_PI);
	  for (int m=0, moffset=0; m<=l; m++) {
	    if (m==0) {
	      in.read((char *)&coefs[loffset+moffset+0][ir], sizeof(double));
	      moffset += 1;
	    } else {
	      in.read((char *)&coefs[loffset+moffset+0][ir], sizeof(double));
	      in.read((char *)&coefs[loffset+moffset+1][ir], sizeof(double));
	      moffset += 2;
	    }
	  }
	}
      }

    } else {

      // Rewind file
      //
      in.clear();
      in.seekg(curpos);

      in.read((char *)&header, sizeof(SphCoefHeader));

      normed = false;

      coefs.resize((header.Lmax+1)*(header.Lmax+1));
      for (auto & v : coefs) v.resize(header.nmax);

      // Read the coefficients
      //
      for (int ir=0; ir<header.nmax; ir++) {
	for (int l=0; l<(header.Lmax+1)*(header.Lmax+1); l++)
	  in.read((char *)&coefs[l][ir], sizeof(double));
      }
    }
  }
  catch (std::istream::failure e) {
    if (not in.eof())
      std::cerr << "Exception reading coefficient file:"
		<< e.what() << std::endl;
    return false;
  }

  // New style norm
  //
  if (not normed) {
    int k = 0;
    for (int l=0; l<=header.Lmax; l++) {
      for (int m=0; m<=l; m++) {
	double fac = sqrt( (0.5*l+0.25)/M_PI *
			   exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );

	if (m != 0) fac *= M_SQRT2;

	// Cosine terms
	for (int ir=0; ir<header.nmax; ir++) coefs[k][ir] *= fac;
	k++;

	// Sine terms
	if (m != 0) {
	  for (int ir=0; ir<header.nmax; ir++) coefs[k][ir] *= fac;
	  k++;
	}
      }
    }

    normed = true;
  }

  coeftomap();

  return true;
}


void SphCoefs::coeftomap()
{
  int k = 0;
  for (int l=0; l<=header.Lmax; l++) {
    for (int m=0; m<=l; m++) {

      LMkey key(l, m);

      cos_c[key].resize(header.nmax);
      if (m != 0) sin_c[key].resize(header.nmax);

      int ik = 0;

      // Cosine terms
      for (int ir=0; ir<header.nmax; ir++) cos_c[key][ir] = coefs[k][ir];
      k++;

      // Sine terms
      if (m != 0) {
	for (int ir=0; ir<header.nmax; ir++) sin_c[key][ir] = coefs[k][ir];
	k++;
      }

    }
  }
}

void SphCoefs::maptocoef()
{
  coefs.resize((header.Lmax+1)*(header.Lmax+1));
  for (auto & v : coefs) v.resize(header.nmax);

  int k = 0;
  for (int l=0; l<=header.Lmax; l++) {
    for (int m=0; m<=l; m++) {

      LMkey key(l, m);

      // Cosine terms
      for (int ir=0; ir<header.nmax; ir++) coefs[k][ir] = cos_c[key][ir];
      k++;

      // Sine terms
      if (m != 0) {
	for (int ir=0; ir<header.nmax; ir++) coefs[k][ir] = sin_c[key][ir];
	k++;
      }
    }
  }
}


void SphCoefs::write(std::ostream& out)
{
  // Convert from map to linear packing for halo coefficient file
  maptocoef();

  // Coefficient magic number
  const unsigned int cmagic = 0xc0a57a2;

  YAML::Node node;

  node["id"    ] = header.id;
  node["time"  ] = header.tnow;
  node["scale" ] = header.scale;
  node["nmax"  ] = header.nmax;
  node["lmax"  ] = header.Lmax;
  node["normed"] = normed;

  // Serialize the node
  //
  YAML::Emitter y; y << node;

  // Get the size of the string
  //
  unsigned int hsize = strlen(y.c_str());

  // Write magic #
  //
  out.write(reinterpret_cast<const char *>(&cmagic),   sizeof(unsigned int));

  // Write YAML string size
  //
  out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));

  // Write YAML string
  //
  out.write(reinterpret_cast<const char *>(y.c_str()), hsize);

  double z;

  // Write the coefficients
  //
  for (int ir=0; ir<header.nmax; ir++) {
    for (int l=0, loffset=0; l<=header.Lmax; loffset+=(2*l+1), l++) {
      for (int m=0, moffset=0; m<=l; m++) {
	if (m==0) {
	  out.write((char *)&(z=coefs[loffset+moffset+0][ir]), sizeof(double));
	  moffset += 1;
	} else {
	  out.write((char *)&(z=coefs[loffset+moffset+0][ir]), sizeof(double));
	  out.write((char *)&(z=coefs[loffset+moffset+1][ir]), sizeof(double));
	  moffset += 2;
	}
      }
    }
  }
}
