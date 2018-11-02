// -*- C++ -*-

#include <Component.H>
#include <Ion.H>
#include <CollideIon.H>
#include <TreeDSMC.H>
#include <EXPException.H>

#include <curand.h>
#include <curand_kernel.h>

#include <cudaUtil.cuH>

// Enable for debugging only
// #define MASS_DEBUG

__constant__ cuFP_t           cuH_H, cuHe_H, cuPH_H, cuPHe_H;
__constant__ cuFP_t           cuH_Emin, cuHe_Emin, cuPH_Emin, cuPHe_Emin;

// Atomic radii in picometers from Clementi, E.; Raimond, D. L.;
// Reinhardt, W. P. (1967). "Atomic Screening Constants from SCF
// Functions. II. Atoms with 37 to 86 Electrons". Journal of Chemical
// Physics 47 (4): 1300-1307.  See also Paper 1, ref. therein.
//
const int numRadii = 87;
__constant__ int cudaRadii[numRadii];

// For construction of evenly spaced interpolation arrays
//
thrust::host_vector<cuFP_t>
resampleArray(const std::vector<cuFP_t>& x, const std::vector<cuFP_t>& y,
	      cuFP_t& dx)
{
  // Get minimum grid spacing
  cuFP_t minH = std::numeric_limits<cuFP_t>::max();
  for (int i=0; i<x.size()-1; i++)
    minH = std::min<cuFP_t>(minH, x[i+1]- x[i]);

  // Resample based on minimum spacing
  int numH = int( (x.back() - x.front())/minH ) + 1;

  thrust::host_vector<cuFP_t> Y(numH);
  
  dx = (x.back() - x.front())/(numH - 1);

  for (int i=0; i<numH; i++) {
    cuFP_t xx = x.front() + dx*i, yy;
    if (xx <= x.front()) {
      yy = y.front();
    } else if (xx >= x.back()) {
      yy = y.back();
    } else {
      auto lb = std::lower_bound(x.begin(), x.end(), xx);
      auto ub = lb;
      if (lb!=x.begin()) lb--;
      if (ub == x.end()) ub = lb--;
      auto a = (*ub - xx)/(*ub - *lb);
      auto b = (xx - *lb)/(*ub - *lb);
      yy = a*y[lb - x.begin()] + b*y[ub - x.begin()];
    }
    Y[i] = yy;
  }

  return Y;
}

// Initialize cross-section look up and interpolation arrays.  Data
// input could be generalized here . . . for later.
//
void CollideIon::cudaElasticInit()
{
  std::vector<int> radii(numRadii, 0);

  radii[1]  =  53;
  radii[2]  =  31;
  radii[3]  =  167;
  radii[4]  =  112;
  radii[5]  =  87;
  radii[6]  =  67;
  radii[7]  =  56;
  radii[8]  =  48;
  radii[9]  =  42;
  radii[10] =  38;
  radii[11] =  190;
  radii[12] =  145;
  radii[13] =  118;
  radii[14] =  111;
  radii[15] =  98;
  radii[16] =  180;
  radii[17] =  79;
  radii[18] =  188;
  radii[19] =  243;
  radii[20] =  194;
  radii[21] =  184;
  radii[22] =  176;
  radii[23] =  171;
  radii[24] =  166;
  radii[25] =  161;
  radii[26] =  156;
  radii[27] =  152;
  radii[28] =  149;
  radii[29] =  145;
  radii[30] =  152;
  radii[31] =  136;
  radii[32] =  125;
  radii[33] =  114;
  radii[34] =  103;
  radii[35] =  94;
  radii[36] =  88;
  radii[37] =  265;
  radii[38] =  219;
  radii[39] =  212;
  radii[40] =  206;
  radii[41] =  198;
  radii[42] =  190;
  radii[43] =  183;
  radii[44] =  178;
  radii[45] =  173;
  radii[46] =  169;
  radii[47] =  172;
  radii[48] =  161;
  radii[49] =  193;
  radii[50] =  217;
  radii[51] =  133;
  radii[52] =  123;
  radii[53] =  198;
  radii[54] =  108;
  radii[55] =  298;
  radii[56] =  268;
  radii[59] =  247;
  radii[60] =  206;
  radii[61] =  205;
  radii[62] =  238;
  radii[63] =  231;
  radii[64] =  233;
  radii[65] =  225;
  radii[66] =  228;
  radii[68] =  226;
  radii[69] =  222;
  radii[70] =  222;
  radii[71] =  217;
  radii[72] =  208;
  radii[73] =  200;
  radii[74] =  193;
  radii[75] =  188;
  radii[76] =  185;
  radii[77] =  180;
  radii[78] =  177;
  radii[79] =  166;
  radii[80] =  171;
  radii[81] =  156;
  radii[82] =  202;
  radii[83] =  143;
  radii[84] =  135;
  radii[86] =  120;

  cuda_safe_call(cudaMemcpyToSymbol(cudaRadii, &radii[0], sizeof(int)*numRadii), 
		 __FILE__, __LINE__, "Error copying cudaRadii");

  // Total cross section from Malik & Trefftz, 1960, Zeitschrift fur Astrophysik, 50, 96-109

  // Column 1 in eV, Column2 in Bohr cross section (pi*a_0^2) units

  std::vector<cuFP_t> eV_H, xs_H;

  eV_H.push_back(0.66356077923727);	xs_H.push_back(28.864);
  eV_H.push_back(0.66576762346358);	xs_H.push_back(29.5088);
  eV_H.push_back(0.71282701193426);	xs_H.push_back(28.2574);
  eV_H.push_back(0.73590227585419);	xs_H.push_back(27.4989);
  eV_H.push_back(0.75936666273882);	xs_H.push_back(26.8542);
  eV_H.push_back(0.80889140369906);	xs_H.push_back(26.3234);
  eV_H.push_back(0.83209592175062);	xs_H.push_back(25.6028);
  eV_H.push_back(0.90677351672548);	xs_H.push_back(24.9204);
  eV_H.push_back(0.95590913472103);	xs_H.push_back(24.2757);
  eV_H.push_back(1.031106467362);	xs_H.push_back(23.745);
  eV_H.push_back(1.05457085424663);	xs_H.push_back(23.1003);
  eV_H.push_back(1.10409559520687);	xs_H.push_back(22.5694);
  eV_H.push_back(1.17903305901478);	xs_H.push_back(21.9628);
  eV_H.push_back(1.22881766880808);	xs_H.push_back(21.5079);
  eV_H.push_back(1.2782131556367);	xs_H.push_back(20.9391);
  eV_H.push_back(1.35379961124236);	xs_H.push_back(20.5221);
  eV_H.push_back(1.45506137966837);	xs_H.push_back(20.1052);
  eV_H.push_back(1.58185287994542);	xs_H.push_back(19.6504);
  eV_H.push_back(1.75999228472356);	xs_H.push_back(19.1958);
  eV_H.push_back(1.91233528596856);	xs_H.push_back(18.7032);
  eV_H.push_back(2.06519530374007);	xs_H.push_back(18.3623);
  eV_H.push_back(2.24360682247953);	xs_H.push_back(17.9835);
  eV_H.push_back(2.42186867854026);	xs_H.push_back(17.5668);
  eV_H.push_back(2.60026659158166);	xs_H.push_back(17.188);
  eV_H.push_back(2.77893661858437);	xs_H.push_back(16.8851);
  eV_H.push_back(2.9830084838763);	xs_H.push_back(16.5064);
  eV_H.push_back(3.21287675270137);	xs_H.push_back(16.1657);
  eV_H.push_back(3.39141072272342);	xs_H.push_back(15.8249);
  eV_H.push_back(3.64683049251644);	xs_H.push_back(15.4463);
  eV_H.push_back(3.87695726960477);	xs_H.push_back(15.1815);
  eV_H.push_back(4.0810291348967);	xs_H.push_back(14.8028);
  eV_H.push_back(4.31091100941984);	xs_H.push_back(14.4621);
  eV_H.push_back(4.54091533522557);	xs_H.push_back(14.1593);
  eV_H.push_back(4.71957175653021);	xs_H.push_back(13.8564);
  eV_H.push_back(4.97525003458648);	xs_H.push_back(13.5537);
  eV_H.push_back(5.28200410318252);	xs_H.push_back(13.1753);
  eV_H.push_back(5.53768238123879);	xs_H.push_back(12.8726);
  eV_H.push_back(5.74227126305723);	xs_H.push_back(12.6456);
  eV_H.push_back(5.97267015410688);	xs_H.push_back(12.4566);
  eV_H.push_back(6.15132657541152);	xs_H.push_back(12.1537);
  eV_H.push_back(6.40726336173105);	xs_H.push_back(11.9268);
  eV_H.push_back(6.61198830053015);	xs_H.push_back(11.7378);
  eV_H.push_back(6.81683569061185);	xs_H.push_back(11.5866);
  eV_H.push_back(6.99562816889715);	xs_H.push_back(11.3216);
  eV_H.push_back(7.20035310769626);	xs_H.push_back(11.1326);
  eV_H.push_back(7.43061594176524);	xs_H.push_back(10.9057);
  eV_H.push_back(7.71209062335465);	xs_H.push_back(10.641);
  eV_H.push_back(7.96789135269351);	xs_H.push_back(10.3762);
  eV_H.push_back(8.27517604351412);	xs_H.push_back(10.1495);
  eV_H.push_back(8.50530282060245);	xs_H.push_back(9.88464);
  eV_H.push_back(8.76123960692198);	xs_H.push_back(9.6578);
  eV_H.push_back(9.06852429774258);	xs_H.push_back(9.43109);
  eV_H.push_back(9.35012143061459);	xs_H.push_back(9.20432);
  eV_H.push_back(9.55484636941369);	xs_H.push_back(9.01526);
  eV_H.push_back(9.78536771174593);	xs_H.push_back(8.86419);
  eV_H.push_back(10.0157666027956);	xs_H.push_back(8.6752);
  eV_H.push_back(10.2717033891151);	xs_H.push_back(8.44835);
  eV_H.push_back(10.5533005219871);	xs_H.push_back(8.22158);
  eV_H.push_back(10.8349112605572);	xs_H.push_back(7.9948);
  eV_H.push_back(11.1421823456797);	xs_H.push_back(7.7681);
  eV_H.push_back(11.4237930842498);	xs_H.push_back(7.54133);
  eV_H.push_back(11.6798523218519);	xs_H.push_back(7.35241);
  eV_H.push_back(11.9615991174026);	xs_H.push_back(7.16356);
  eV_H.push_back(12.2176583550047);	xs_H.push_back(6.97464);
  eV_H.push_back(12.4223832938038);	xs_H.push_back(6.78558);
  eV_H.push_back(12.7041164836565);	xs_H.push_back(6.59673);
  eV_H.push_back(12.9858496735092);	xs_H.push_back(6.40788);
  eV_H.push_back(13.2163710158414);	xs_H.push_back(6.25682);
  eV_H.push_back(13.4212320116212);	xs_H.push_back(6.10568);
  eV_H.push_back(13.600541506433);	xs_H.push_back(5.9924);;

  cuFP_t dx;

  xsc_H = resampleArray(eV_H, xs_H, dx);

  cuda_safe_call(cudaMemcpyToSymbol(cuH_Emin, &eV_H[0], sizeof(cuFP_t)),
		 __FILE__, __LINE__, "Error copying cuH_Emin");

  cuda_safe_call(cudaMemcpyToSymbol(cuH_H, &dx, sizeof(cuFP_t)),
		 __FILE__, __LINE__, "Error copying cuH_H");

  // Total cross section from LaBahn & Callaway, 1966, Phys. Rev., 147, 50, 28-40
  //
  std::vector<cuFP_t> eV_He, xs_He;

  eV_He.push_back(0.0972135);	xs_He.push_back(62.2773619684373);
  eV_He.push_back(0.212908);	xs_He.push_back(65.3193661349083);
  eV_He.push_back(0.251768);	xs_He.push_back(66.6419766420696);
  eV_He.push_back(0.440878);	xs_He.push_back(67.8332685763108);
  eV_He.push_back(0.704798);	xs_He.push_back(68.6287198361997);
  eV_He.push_back(1.11846);	xs_He.push_back(68.7641224795695);
  eV_He.push_back(1.5694);	xs_He.push_back(68.5696578943123);
  eV_He.push_back(1.86971);	xs_He.push_back(68.109100411296 );
  eV_He.push_back(2.20762);	xs_He.push_back(67.6491712468105);
  eV_He.push_back(2.50774);	xs_He.push_back(66.9903792673527);
  eV_He.push_back(2.77027);	xs_He.push_back(66.3312731286295);
  eV_He.push_back(3.07045);	xs_He.push_back(65.7387687541625);
  eV_He.push_back(3.25779);	xs_He.push_back(65.0790342969086);
  eV_He.push_back(3.44532);	xs_He.push_back(64.6178484953617);
  eV_He.push_back(3.78297);	xs_He.push_back(63.8930830701785);
  eV_He.push_back(4.0079);	xs_He.push_back(63.2339769314554);
  eV_He.push_back(4.2329);	xs_He.push_back(62.6405300791922);
  eV_He.push_back(4.57055);	xs_He.push_back(61.9160788132744);
  eV_He.push_back(4.79555);	xs_He.push_back(61.3229461202767);
  eV_He.push_back(5.02048);	xs_He.push_back(60.6635258222882);
  eV_He.push_back(5.20782);	xs_He.push_back(60.0041055242997);
  eV_He.push_back(5.50788);	xs_He.push_back(59.2790259398512);
  eV_He.push_back(5.883);	xs_He.push_back(58.4226277824826);
  eV_He.push_back(6.14552);	xs_He.push_back(57.7635216437595);
  eV_He.push_back(6.48318);	xs_He.push_back(57.0390703778416);
  eV_He.push_back(6.89576);	xs_He.push_back(56.0507253290223);
  eV_He.push_back(7.27088);	xs_He.push_back(55.1943271716537);
  eV_He.push_back(7.68347);	xs_He.push_back(54.2059821228344);
  eV_He.push_back(7.98353);	xs_He.push_back(53.4809025383858);
  eV_He.push_back(8.32118);	xs_He.push_back(52.756451272468);
  eV_He.push_back(8.6963);	xs_He.push_back(51.9000531150995);
  eV_He.push_back(8.99617);	xs_He.push_back(50.9767390342094);
  eV_He.push_back(9.33408);	xs_He.push_back(50.5168098697239);
  eV_He.push_back(9.67173);	xs_He.push_back(49.7920444445407);
  eV_He.push_back(9.97191);	xs_He.push_back(49.1995400700737);
  eV_He.push_back(10.2346);	xs_He.push_back(48.6726949820667);
  eV_He.push_back(10.4596);	xs_He.push_back(48.0795622890689);
  eV_He.push_back(10.7222);	xs_He.push_back(47.5527172010619);
  eV_He.push_back(10.9849);	xs_He.push_back(47.0921597180456);
  eV_He.push_back(11.2852);	xs_He.push_back(46.565628789304 );
  eV_He.push_back(11.5478);	xs_He.push_back(46.038783701297 );
  eV_He.push_back(11.773);	xs_He.push_back(45.57759789975  );
  eV_He.push_back(12.0731);	xs_He.push_back(44.9191200795576);
  eV_He.push_back(12.3734);	xs_He.push_back(44.4585625965413);
  eV_He.push_back(12.7488);	xs_He.push_back(43.866686540605 );
  eV_He.push_back(13.0489);	xs_He.push_back(43.2738680068726);
  eV_He.push_back(13.3867);	xs_He.push_back(42.6153901866802);
  eV_He.push_back(13.7996);	xs_He.push_back(41.9578548442838);
  eV_He.push_back(14.1375);	xs_He.push_back(41.4976115205329);
  eV_He.push_back(14.4752);	xs_He.push_back(40.9054213053313);
  eV_He.push_back(14.8132);	xs_He.push_back(40.5114655865711);
  eV_He.push_back(15.1509);	xs_He.push_back(39.8529877663787);
  eV_He.push_back(15.4889);	xs_He.push_back(39.4590320476185);
  eV_He.push_back(15.7893);	xs_He.push_back(39.064762169593 );
  eV_He.push_back(16.1271);	xs_He.push_back(38.5385454001167);
  eV_He.push_back(16.465);	xs_He.push_back(38.0123286306404);
  eV_He.push_back(16.8404);	xs_He.push_back(37.4864260204295);
  eV_He.push_back(17.1407);	xs_He.push_back(37.0258685374132);
  eV_He.push_back(17.5162);	xs_He.push_back(36.566253532193 );
  eV_He.push_back(17.8917);	xs_He.push_back(36.1063243677075);
  eV_He.push_back(18.2672);	xs_He.push_back(35.6463952032219);
  eV_He.push_back(18.6426);	xs_He.push_back(35.120492593011 );
  eV_He.push_back(19.0181);	xs_He.push_back(34.6608775877908);
  eV_He.push_back(19.3561);	xs_He.push_back(34.2669218690307);
  eV_He.push_back(19.6941);	xs_He.push_back(33.8729661502705);
  eV_He.push_back(20.0321);	xs_He.push_back(33.478696272245);
  eV_He.push_back(20.3324);	xs_He.push_back(33.0844263942195);
  eV_He.push_back(20.708);	xs_He.push_back(32.6907848347247);
  eV_He.push_back(21.046);	xs_He.push_back(32.2968291159645);
  eV_He.push_back(21.3839);	xs_He.push_back(31.836899951479 );
  eV_He.push_back(21.7594);	xs_He.push_back(31.37700220292  );
  eV_He.push_back(22.1725);	xs_He.push_back(30.9174814454794);
  eV_He.push_back(22.548);	xs_He.push_back(30.523808470058 );
  eV_He.push_back(22.9612);	xs_He.push_back(30.1304182379755);
  eV_He.push_back(23.3369);	xs_He.push_back(29.8689434814172);
  eV_He.push_back(23.5995);	xs_He.push_back(29.3421612252633);
  eV_He.push_back(23.9752);	xs_He.push_back(29.080686468705 );
  eV_He.push_back(24.3506);	xs_He.push_back(28.4886847490626);
  eV_He.push_back(24.8389);	xs_He.push_back(28.0958914195842);
  eV_He.push_back(25.2144);	xs_He.push_back(27.6360879188048);
  eV_He.push_back(25.8906);	xs_He.push_back(27.1125729190106);
  eV_He.push_back(26.3788);	xs_He.push_back(26.5875499547427);
  eV_He.push_back(26.7543);	xs_He.push_back(26.1277464539633);
  eV_He.push_back(27.2427);	xs_He.push_back(25.734953124485 );
  eV_He.push_back(27.7686);	xs_He.push_back(25.342473954272 );
  eV_He.push_back(28.2944);	xs_He.push_back(24.8177337333429);
  eV_He.push_back(28.8205);	xs_He.push_back(24.5574841979195);
  eV_He.push_back(29.2337);	xs_He.push_back(24.164093965837 );
  eV_He.push_back(29.8723);	xs_He.push_back(23.706363916209 );
  eV_He.push_back(30.3607);	xs_He.push_back(23.3135705867306);
  eV_He.push_back(30.8868);	xs_He.push_back(23.1194201607388);
  eV_He.push_back(31.4127);	xs_He.push_back(22.7269409905258);
  eV_He.push_back(31.9388);	xs_He.push_back(22.4666600391759);
  eV_He.push_back(32.615);	xs_He.push_back(21.9431450393817);
  eV_He.push_back(33.3665);	xs_He.push_back(21.5524251610547);
  eV_He.push_back(34.0429);	xs_He.push_back(21.2272389054817);
  eV_He.push_back(34.6064);	xs_He.push_back(20.8350424786075);
  eV_He.push_back(35.0948);	xs_He.push_back(20.5083796744872);
  eV_He.push_back(35.5457);	xs_He.push_back(20.2475018205331);
  eV_He.push_back(36.0342);	xs_He.push_back(20.0530372352759);
  eV_He.push_back(36.4851);	xs_He.push_back(19.7921907972484);
  eV_He.push_back(37.0112);	xs_He.push_back(19.59800895533  );
  eV_He.push_back(37.5372);	xs_He.push_back(19.2716288945485);
  eV_He.push_back(38.0258);	xs_He.push_back(19.0771957252179);
  eV_He.push_back(38.5519);	xs_He.push_back(18.8830138832995);
  eV_He.push_back(38.9277);	xs_He.push_back(18.6876696520993);
  eV_He.push_back(39.4537);	xs_He.push_back(18.4273887007494);
  eV_He.push_back(39.9423);	xs_He.push_back(18.2329555314187);
  eV_He.push_back(40.506);	xs_He.push_back(18.0390878487657);
  eV_He.push_back(40.9569);	xs_He.push_back(17.7782099948116);
  eV_He.push_back(41.483);	xs_He.push_back(17.5840595688197);
  eV_He.push_back(41.9715);	xs_He.push_back(17.3895949835625);
  eV_He.push_back(42.3472);	xs_He.push_back(17.1281516429308);
  eV_He.push_back(42.7606);	xs_He.push_back(16.9991892645009);
  eV_He.push_back(43.174);	xs_He.push_back(16.8702583019976);
  eV_He.push_back(43.5498);	xs_He.push_back(16.6749140707974);
  eV_He.push_back(43.9631);	xs_He.push_back(16.479852582936);
  eV_He.push_back(44.339);	xs_He.push_back(16.3506074611673);
  eV_He.push_back(44.6772);	xs_He.push_back(16.2210795960598);
  eV_He.push_back(45.2033);	xs_He.push_back(16.0268977541414);
  eV_He.push_back(45.504);	xs_He.push_back(15.9631862551266);
  eV_He.push_back(45.8798);	xs_He.push_back(15.8339411333579);
  eV_He.push_back(46.2556);	xs_He.push_back(15.7046960115892);
  eV_He.push_back(46.5939);	xs_He.push_back(15.5751681464817);
  eV_He.push_back(46.857);	xs_He.push_back(15.5111424882016);
  eV_He.push_back(47.1952);	xs_He.push_back(15.3816146230941);
  eV_He.push_back(47.571);	xs_He.push_back(15.2523695013254);
  eV_He.push_back(47.8717);	xs_He.push_back(15.188626586384 );
  eV_He.push_back(48.1348);	xs_He.push_back(15.1246009281039);
  eV_He.push_back(48.4354);	xs_He.push_back(14.9947903196575);
  eV_He.push_back(48.7736);	xs_He.push_back(14.8652310386235);
  eV_He.push_back(49.0368);	xs_He.push_back(14.8673359057014);
  eV_He.push_back(49.2622);	xs_He.push_back(14.7368969787244);
  eV_He.push_back(49.5254);	xs_He.push_back(14.7389704298757);
  eV_He.push_back(49.7885);	xs_He.push_back(14.6749447715956);
  eV_He.push_back(49.9763);	xs_He.push_back(14.4781239918482);

  xsc_He = resampleArray(eV_He, xs_He, dx);

  cuda_safe_call(cudaMemcpyToSymbol(cuH_Emin, &eV_He[0], sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuHe_Emin");

  cuda_safe_call(cudaMemcpyToSymbol(cuH_H, &dx, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuHe_H");

  // Interpolated from Figure 1 of "Elastic scattering and charge
  // transfer in slow collisions: isotopes of H and H + colliding with
  // isotopes of H and with He" by Predrag S Krstić and David R Schultz,
  // 1999 J. Phys. B: At. Mol. Opt. Phys. 32 3485
  //

  std::vector<cuFP_t> eV_pH, xs_pH;

  eV_pH.push_back(-0.994302);	xs_pH.push_back(2.86205);
  eV_pH.push_back(-0.897482);	xs_pH.push_back(2.90929);
  eV_pH.push_back(-0.801179);	xs_pH.push_back(2.86016);
  eV_pH.push_back(-0.691555);	xs_pH.push_back(2.89417);
  eV_pH.push_back(-0.588753);	xs_pH.push_back(2.85638);
  eV_pH.push_back(-0.49242);	xs_pH.push_back(2.81291);
  eV_pH.push_back(-0.395965);	xs_pH.push_back(2.79213);
  eV_pH.push_back(-0.292839);	xs_pH.push_back(2.8148);
  eV_pH.push_back(-0.19019);	xs_pH.push_back(2.74866);
  eV_pH.push_back(-0.0872765);	xs_pH.push_back(2.73165);
  eV_pH.push_back(0.00935082);	xs_pH.push_back(2.74299);
  eV_pH.push_back(0.112152);	xs_pH.push_back(2.7052);
  eV_pH.push_back(0.208688);	xs_pH.push_back(2.69953);
  eV_pH.push_back(0.311612);	xs_pH.push_back(2.68441);
  eV_pH.push_back(0.401578);	xs_pH.push_back(2.65417);
  eV_pH.push_back(0.517468);	xs_pH.push_back(2.65606);
  eV_pH.push_back(0.613862);	xs_pH.push_back(2.62394);
  eV_pH.push_back(0.716846);	xs_pH.push_back(2.62016);
  eV_pH.push_back(0.819688);	xs_pH.push_back(2.58992);
  eV_pH.push_back(0.909797);	xs_pH.push_back(2.58614);
  eV_pH.push_back(1.01906);	xs_pH.push_back(2.55213);
  eV_pH.push_back(1.1092);	xs_pH.push_back(2.55402);
  eV_pH.push_back(1.21203);	xs_pH.push_back(2.52189);
  eV_pH.push_back(1.3085);	xs_pH.push_back(2.50488);
  eV_pH.push_back(1.41149);	xs_pH.push_back(2.5011);
  eV_pH.push_back(1.52077);	xs_pH.push_back(2.47087);
  eV_pH.push_back(1.61715);	xs_pH.push_back(2.43685);
  eV_pH.push_back(1.71368);	xs_pH.push_back(2.42929);
  eV_pH.push_back(1.81666);	xs_pH.push_back(2.42551);
  eV_pH.push_back(1.9131);	xs_pH.push_back(2.40094);
  eV_pH.push_back(2.0159);	xs_pH.push_back(2.36315);

  xsc_pH = resampleArray(eV_pH, xs_pH, dx);

  cuda_safe_call(cudaMemcpyToSymbol(cuPH_Emin, &eV_pH[0], sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuPH_Emin");

  cuda_safe_call(cudaMemcpyToSymbol(cuPH_H, &dx, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuPH_H");

  // Interpolated from the top panel of Figure 4, op. cit.
  //
  std::vector<cuFP_t> eV_pHe, xs_pHe;

  eV_pHe.push_back(-0.984127);	xs_pHe.push_back(2.68889);
  eV_pHe.push_back(-0.904762);	xs_pHe.push_back(2.68889);
  eV_pHe.push_back(-0.825397);	xs_pHe.push_back(2.68889);
  eV_pHe.push_back(-0.753968);	xs_pHe.push_back(2.64444);
  eV_pHe.push_back(-0.674603);	xs_pHe.push_back(2.6);
  eV_pHe.push_back(-0.595238);	xs_pHe.push_back(2.57778);
  eV_pHe.push_back(-0.515873);	xs_pHe.push_back(2.57778);
  eV_pHe.push_back(-0.444444);	xs_pHe.push_back(2.55556);
  eV_pHe.push_back(-0.373016);	xs_pHe.push_back(2.48889);
  eV_pHe.push_back(-0.293651);	xs_pHe.push_back(2.44444);
  eV_pHe.push_back(-0.214286);	xs_pHe.push_back(2.46667);
  eV_pHe.push_back(-0.142857);	xs_pHe.push_back(2.44444);
  eV_pHe.push_back(-0.0634921);	xs_pHe.push_back(2.4);
  eV_pHe.push_back(0.015873);	xs_pHe.push_back(2.37778);
  eV_pHe.push_back(0.0952381);	xs_pHe.push_back(2.37778);
  eV_pHe.push_back(0.166667);	xs_pHe.push_back(2.33333);
  eV_pHe.push_back(0.246032);	xs_pHe.push_back(2.28889);
  eV_pHe.push_back(0.325397);	xs_pHe.push_back(2.28889);
  eV_pHe.push_back(0.404762);	xs_pHe.push_back(2.28889);
  eV_pHe.push_back(0.47619);	xs_pHe.push_back(2.24444);
  eV_pHe.push_back(0.555556);	xs_pHe.push_back(2.2);
  eV_pHe.push_back(0.634921);	xs_pHe.push_back(2.17778);
  eV_pHe.push_back(0.706349);	xs_pHe.push_back(2.2);
  eV_pHe.push_back(0.785714);	xs_pHe.push_back(2.17778);
  eV_pHe.push_back(0.865079);	xs_pHe.push_back(2.13333);
  eV_pHe.push_back(0.936508);	xs_pHe.push_back(2.08889);
  eV_pHe.push_back(1.01587);	xs_pHe.push_back(2.06667);
  eV_pHe.push_back(1.09524);	xs_pHe.push_back(2.08889);
  eV_pHe.push_back(1.16667);	xs_pHe.push_back(2.06667);
  eV_pHe.push_back(1.24603);	xs_pHe.push_back(2.04444);
  eV_pHe.push_back(1.3254);	xs_pHe.push_back(2.02222);
  eV_pHe.push_back(1.40476);	xs_pHe.push_back(1.97778);
  eV_pHe.push_back(1.47619);	xs_pHe.push_back(1.93333);
  eV_pHe.push_back(1.55556);	xs_pHe.push_back(1.91111);
  eV_pHe.push_back(1.63492);	xs_pHe.push_back(1.91111);
  eV_pHe.push_back(1.71429);	xs_pHe.push_back(1.91111);
  eV_pHe.push_back(1.79365);	xs_pHe.push_back(1.91111);
  eV_pHe.push_back(1.87302);	xs_pHe.push_back(1.91111);
  eV_pHe.push_back(1.95238);	xs_pHe.push_back(1.91111);

  xsc_pHe = resampleArray(eV_pHe, xs_pHe, dx);

  cuda_safe_call(cudaMemcpyToSymbol(cuPHe_Emin, &eV_pHe[0], sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuPHe_Emin");

  cuda_safe_call(cudaMemcpyToSymbol(cuPHe_H, &dx, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuPHe_H");
}

__device__
cuFP_t cudaGeometric(int Z)
{
  if (Z>0 and Z< numRadii) {
    return cudaRadii[Z] * 1.0e-3;
  } else {
    return 0.0;
  }
}
		 
__device__
cuFP_t cudaElasticInterp(cuFP_t E, cuFP_t Emin, cuFP_t H, dArray<cuFP_t> xsc,
			 bool pin = true)
{
  int indx = 0;
  if (E >= Emin+H*xsc._s) indx = xsc._s - 2;
  if (E <  Emin)          indx = 0;
  else                    indx = floor( (E - Emin)/H );

  cuFP_t a = (E - Emin - H*(indx+0))/H;
  cuFP_t b = (Emin + H*(indx+1) - E)/H;

  // Enforce return value to grid boundaries for off-grid ordinates.
  // Otherwise, values will be extrapolated.
  if (pin) {
    if (a < 0.0) return xsc._v[0];
    if (b < 0.0) return xsc._v[xsc._s-1];
  }

  return a*xsc._v[indx] + b*xsc._v[indx+1];
}


// Global symbols for coordinate transformation
//
__device__ __constant__
cuFP_t ionEminGrid, ionEmaxGrid, ionDeltaEGrid;

__device__ __constant__
int ionEgridNumber, ionRadRecombNumber;

__global__
void testConstantsIon()
{
  printf("** Egrid(min) = %f\n", ionEminGrid);
  printf("** Egrid(max) = %f\n", ionEmaxGrid);
  printf("** Egrid(del) = %f\n", ionDeltaEGrid);
  printf("** Egrid(num) = %d\n", ionEgridNumber);
  printf("** Rgrid(num) = %d\n", ionRadRecombNumber);
}

void chdata::cuda_initialize_textures()
{
  size_t ionSize = IonList.size();

  // Interpolation data array
  //
  cuF0array.resize(ionSize, 0);
  cuFFarray.resize(ionSize, 0);
  cuRCarray.resize(ionSize, 0);
  cuCEarray.resize(ionSize, 0);
  cuCIarray.resize(ionSize, 0);
  cuPIarray.resize(ionSize   );

  // Texture object array
  //
  cuIonElem.resize(ionSize);

  // Total photo-ionization rate
  //
  std::vector<cuFP_t> phRate(ionSize, 0.0);

  size_t k = 0;

  for (auto v : IonList) {

    IonPtr I = v.second;
    cuIonElement& E = cuIonElem[k];

    // The free-free array
    if (E.C>1) {
      cudaTextureDesc texDesc;

      memset(&texDesc, 0, sizeof(texDesc));
      texDesc.readMode = cudaReadModeElementType;
      texDesc.filterMode = cudaFilterModePoint;
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.addressMode[1] = cudaAddressModeClamp;
      texDesc.addressMode[2] = cudaAddressModeClamp;
  
      // Temporary storage
      //
      std::vector<cuFP_t> h_buffer0(I->NfreeFreeGrid, 0.0);

      cuFP_t *d_Interp;

      cuda_safe_call(cudaMalloc((void **)&d_Interp, I->NfreeFreeGrid*CHCUMK*sizeof(cuFP_t)),
		     __FILE__, __LINE__,
		     "Error allocating d_Interp1 for texture construction");
  
      std::vector<cuFP_t> h_buffer1(I->NfreeFreeGrid*CHCUMK, 0.0);

      double delC = 1.0/(CHCUMK-1);

      // Copy cross section values to buffer
      //
      for (int i = 0; i < I->NfreeFreeGrid; i++) {

	h_buffer0[i] = I->freeFreeGrid[i].back();
	
	// Unit normalized cumulative distribution
	//
	size_t tsize = I->freeFreeGrid[i].size();
	std::vector<double> temp(tsize);
	for (int j = 0; j < tsize; j++) {	
	  temp[j] = I->freeFreeGrid[i][j]/h_buffer0[i];
	}

	// End points
	//
	h_buffer1[i                              ] = I->kgrid[0];
	h_buffer1[i + (CHCUMK-1)*I->NfreeFreeGrid] = I->kgrid[tsize-1];

	// Remap to even grid
	//
	for (int j=1; j<CHCUMK-1; j++) {

	  double C = delC*j;

	  // Points to first element that is not < C
	  // but may be equal
	  std::vector<double>::iterator lb = 
	    std::lower_bound(temp.begin(), temp.end(), C);
    
	  // Assign upper end of range to the
	  // found element
	  //
	  std::vector<double>::iterator ub = lb;
	  //
	  // If is the first element, increment
	  // the upper boundary
	  //
	  if (lb == temp.begin()) { if (temp.size()>1) ub++; }
	  //
	  // Otherwise, decrement the lower boundary
	  //
	  else { lb--; }
    
	  // Compute the associated indices
	  //
	  size_t ii = lb - temp.begin();
	  size_t jj = ub - temp.begin();
	  double kk = I->kgrid[ii];
	  
	  // Linear interpolation
	  //
	  if (*ub > *lb) {
	    double d = *ub - *lb;
	    double a = (C - *lb) / d;
	    double b = (*ub - C) / d;
	    /*
	    std::cout << "[a, b] = [" << a << ", " << b << "]"
		      << ", c = " << C << std::endl;
	    */
	    kk  = a * I->kgrid[ii] + b * I->kgrid[jj];
	  }

	  h_buffer1[i + j*I->NfreeFreeGrid] = kk;

	} // END: cumululative array loop

      } // END: energy loop

      // Copy 1-dim data to device
      //
      size_t tsize = I->NfreeFreeGrid*sizeof(cuFP_t);

      // Allocate CUDA array in device memory (a one-dimension 'channel')
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<int2>();
#endif
      
      // std::cout << "Allocating cuF0array[" << k << "]" << std::endl;
      cuda_safe_call(cudaMallocArray(&cuF0array[k], &channelDesc1, I->NfreeFreeGrid), __FILE__, __LINE__, "malloc cuArray");

      cuda_safe_call(cudaMemcpyToArray(cuF0array[k], 0, 0, &h_buffer0[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

      // Specify 1-d texture

      cudaResourceDesc resDesc;

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array = cuF0array[k];

      cuda_safe_call(cudaCreateTextureObject(&E.ff_0, &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");

      // Copy data to device
      tsize = I->NfreeFreeGrid*CHCUMK*sizeof(cuFP_t);
      cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer1[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");
    
      // cuda 2d Array Descriptor
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<int2>();
#endif
      // cuda 2d Array
      //
      cuda_safe_call(cudaMalloc3DArray(&cuFFarray[k], &channelDesc2, make_cudaExtent(I->NfreeFreeGrid, CHCUMK, 1), 0), __FILE__, __LINE__, "Error allocating cuArray for 3d texture");
      
      // Array creation
      //
      cudaMemcpy3DParms copyParams = {0};
  
      copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, I->NfreeFreeGrid*sizeof(cuFP_t), I->NfreeFreeGrid, CHCUMK);
      copyParams.dstArray = cuFFarray[k];
      copyParams.extent   = make_cudaExtent(I->NfreeFreeGrid, CHCUMK, 1);
      copyParams.kind     = cudaMemcpyDeviceToDevice;
      
      cuda_safe_call(cudaMemcpy3D(&copyParams), __FILE__, __LINE__, "Error in copying 3d pitched array");

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array  = cuFFarray[k];
    
      cuda_safe_call
	(cudaCreateTextureObject(&E.ff_d, &resDesc, &texDesc, NULL),
	 __FILE__, __LINE__, "Failure in 2d texture creation");
      
      cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");
    }

    // Radiative recombination texture (1-d)
    //
    if (E.C>1) {
      // Allocate CUDA array in device memory (a one-dimension 'channel')
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
    
      // Size of interpolation array
      //
      size_t tsize = I->NradRecombGrid*sizeof(cuFP_t);

      cudaTextureDesc texDesc;
      
      memset(&texDesc, 0, sizeof(cudaTextureDesc));
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.filterMode = cudaFilterModePoint;
      texDesc.readMode = cudaReadModeElementType;
      texDesc.normalizedCoords = 0;
      
      thrust::host_vector<cuFP_t> tt(I->NradRecombGrid);
      
      cuda_safe_call(cudaMallocArray(&cuRCarray[k], &channelDesc, I->NradRecombGrid), __FILE__, __LINE__, "malloc cuArray");

      // Copy to device memory some data located at address h_data
      // in host memory
      for (size_t n = 0; n < I->NradRecombGrid; n++) tt[n] = I->radRecombGrid[n];
    
      cuda_safe_call(cudaMemcpyToArray(cuRCarray[k], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

      // Specify texture

      cudaResourceDesc resDesc;

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array = cuRCarray[k];
      
      cuda_safe_call(cudaCreateTextureObject(&E.rc_d, &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
    }

    // The collisional excitation array

    if (E.C <= E.Z) {

      E.ceEmin = I->collideEmin;
      E.ceEmax = I->collideEmax;
      E.ceDelE = I->delCollideE;
      E.NColl  = I->NcollideGrid;

      /*
      std::cout << " k=" << k
		<< " Emin=" << E.ceEmin
		<< " Emax=" << E.ceEmax
		<< " delE=" << E.ceDelE
		<< std::endl;
      */

      cudaTextureDesc texDesc;
      
      memset(&texDesc, 0, sizeof(texDesc));
      texDesc.readMode = cudaReadModeElementType;
      texDesc.filterMode = cudaFilterModePoint;
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.addressMode[1] = cudaAddressModeClamp;
      texDesc.addressMode[2] = cudaAddressModeClamp;
  
      // Temporary storage
      //
      cuFP_t *d_Interp;
      /*
      std::cout << "Size(" << I->Z << ", " << I->C << ")="
		<< I->NcollideGrid << std::endl;
      */
      cuda_safe_call(cudaMalloc((void **)&d_Interp, I->NcollideGrid*2*sizeof(cuFP_t)),
		     __FILE__, __LINE__,
		     "Error allocating d_Interp for texture construction");
  
      std::vector<cuFP_t> h_buffer(I->NcollideGrid*2, 0.0);

      // Copy vectors to buffer
      //
      for (int i = 0; i < I->NcollideGrid; i++) {
	h_buffer[i                  ] = I->collideDataGrid[i].back().first;
	h_buffer[i + I->NcollideGrid] = I->collideDataGrid[i].back().second;
      }
      
      // Copy data to device
      cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer[0], I->NcollideGrid*2*sizeof(cuFP_t), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");
    
      // cudaArray Descriptor
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
      // cuda Array
      //
      cuda_safe_call(cudaMalloc3DArray(&cuCEarray[k], &channelDesc, make_cudaExtent(I->NcollideGrid, 2, 1), 0), __FILE__, __LINE__, "Error allocating cuArray for 3d texture");
    
      // Array creation
      //
      cudaMemcpy3DParms copyParams = {0};
      
      copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, I->NcollideGrid*sizeof(cuFP_t), I->NcollideGrid, 2);
      copyParams.dstArray = cuCEarray[k];
      copyParams.extent   = make_cudaExtent(I->NcollideGrid, 2, 1);
      copyParams.kind     = cudaMemcpyDeviceToDevice;
      
      cuda_safe_call(cudaMemcpy3D(&copyParams), __FILE__, __LINE__, "Error in copying 3d pitched array");
      
      cudaResourceDesc resDesc;
      
      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array  = cuCEarray[k];
      
      cuda_safe_call
	(cudaCreateTextureObject(&E.ce_d, &resDesc, &texDesc, NULL),
	 __FILE__, __LINE__, "Failure in 2d texture creation");
      
      cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");
    }

    if (E.C <= E.Z) {

      E.ciEmin = I->ionizeEmin;
      E.ciEmax = I->ionizeEmax;
      E.ciDelE = I->DeltaEGrid;
      E.NIonz  = I->NionizeGrid;

      // Allocate CUDA array in device memory (a one-dimension 'channel')
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
    
      // Size of interpolation array
      //
      size_t tsize = I->NionizeGrid*sizeof(cuFP_t);
      
      cudaTextureDesc texDesc;

      memset(&texDesc, 0, sizeof(cudaTextureDesc));
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.filterMode = cudaFilterModePoint;
      texDesc.readMode = cudaReadModeElementType;
      texDesc.normalizedCoords = 0;
      
      thrust::host_vector<cuFP_t> tt(I->NionizeGrid);
      
      /*
      std::cout << "Size(" << I->Z << ", " << I->C << ")="
		<< I->NionizeGrid << std::endl;
      */

      cuda_safe_call(cudaMallocArray(&cuCIarray[k], &channelDesc, I->NionizeGrid), __FILE__, __LINE__, "malloc cuArray");

      // Copy to device memory some data located at address h_data
      // in host memory
      for (size_t n = 0; n < I->NionizeGrid; n++) tt[n] = I->ionizeDataGrid[n];
      
      cuda_safe_call(cudaMemcpyToArray(cuCIarray[k], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");
      
      // Specify texture

      cudaResourceDesc resDesc;

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array = cuCIarray[k];
      
      cuda_safe_call(cudaCreateTextureObject(&E.ci_d, &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");

      // Photoionization array
      //
      if (I->ib_type != Ion::none) {

	thrust::host_vector<cuFP_t> piCum(CHCUMK, 0.0);
	piCum[CHCUMK-1] = 1.0;
      
	double delC = 1.0/(CHCUMK-1);
      
	if (not I->IBinit) I->IBcreate();
      
	E.piTotl = I->IBtotl;

	// Copy cross section values to buffer
	//
	for (int j=1; j<CHCUMK-1; j++) {

	  // Location in cumulative cross section grid
	  //
	  double C = delC*j;

	  // Interpolate the cross section array
	  //
	
	  // Points to first element that is not < rn
	  // but may be equal
	  std::vector<double>::iterator lb = 
	    std::lower_bound(I->IBcum.begin(), I->IBcum.end(), C);
	  
	  // Assign upper end of range to the
	  // found element
	  //
	  std::vector<double>::iterator ub = lb;
	  //
	  // If is the first element, increment
	  // the upper boundary
	  //
	  if (lb == I->IBcum.begin()) { if (I->IBcum.size()>1) ub++; }
	  //
	  // Otherwise, decrement the lower boundary
	  //
	  else { lb--; }
	  
	  // Compute the associated indices
	  //
	  size_t ii = lb - I->IBcum.begin();
	  size_t jj = ub - I->IBcum.begin();
	  double nu = I->nugrid[ii];
	  
	  // Linear interpolation
	  //
	  if (*ub > *lb) {
	    double d = *ub - *lb;
	    double a = (C - *lb) / d;
	    double b = (*ub - C) / d;
	    nu  = a * I->nugrid[ii] + b * I->nugrid[jj];
	  }
	  
	  piCum[j] = (nu - 1.0)*I->ip;
	}
	
	// std::cout << "Allocating pi_0[" << k << "]" << std::endl;

	// Create storage on device
	cuPIarray[k] = piCum;

	// Assign pointer
	E.pi_0 = thrust::raw_pointer_cast(&cuPIarray[k][0]);

      } // END: cumululative array loop

    } // END: ions with electrons
    
    // Increment counter
    k++;	
    
  } // END: IonList

}

void chdata::cuda_initialize_grid_constants()
{
  double Emin, Emax, delE;
  int NE, NR;

  for (auto v : IonList) {
    Emin = v.second->EminGrid;
    Emax = v.second->EmaxGrid;
    delE = v.second->DeltaEGrid;

    NE   = v.second->NfreeFreeGrid;

    if (v.first.second>1) {
      NR = v.second->NradRecombGrid;
      break;
    }
  }

  cuFP_t f;

  // Copy constants to device
  //
  cuda_safe_call(cudaMemcpyToSymbol(ionEminGrid, &(f=Emin),
				    sizeof(cuFP_t), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionEminGrid");

  cuda_safe_call(cudaMemcpyToSymbol(ionEmaxGrid, &(f=Emax),
				    sizeof(cuFP_t), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionEmaxGrid");

  cuda_safe_call(cudaMemcpyToSymbol(ionDeltaEGrid, &(f=delE),
				    sizeof(cuFP_t), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionDeltaEGrid");

  cuda_safe_call(cudaMemcpyToSymbol(ionEgridNumber, &NE,
				    sizeof(int), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionEgridNumber");

  cuda_safe_call(cudaMemcpyToSymbol(ionRadRecombNumber, &NR,
				    sizeof(int), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionRadRecombNumber");
}


__device__
void computeFreeFree
(cuFP_t E, cuFP_t rr, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem)
{
  // value of h-bar * c in eV*nm
  //
  constexpr double hbc = 197.327;

  // Enforce minimum and maximum energies
  //
  if (E<ionEminGrid) E = ionEminGrid;
  if (E>ionEmaxGrid) E = ionEmaxGrid;

  size_t indx = std::floor( (E - ionEminGrid)/ionDeltaEGrid );
    
  if (indx >= ionEgridNumber - 1) indx = ionEgridNumber-2;

  double eA = ionEminGrid + ionDeltaEGrid*indx;
  double eB = ionEminGrid + ionDeltaEGrid*(indx+1);
  
  double A = (eB - E)/ionDeltaEGrid;
  double B = (E - eA)/ionDeltaEGrid;
  
  // Location in cumulative cross section grid
  //
  double rn = rr;
  double dC = 1.0/(CHCUMK-1);
  int lb    = rn/dC;
  cuFP_t k[4];

  // Interpolate the cross section array
  //
#if cuREAL == 4
  k[0]  = tex3D<float>(elem.ff_d, indx,   lb  , 0);
  k[1]  = tex3D<float>(elem.ff_d, indx+1, lb  , 0);
  k[2]  = tex3D<float>(elem.ff_d, indx,   lb+1, 0);
  k[3]  = tex3D<float>(elem.ff_d, indx+1, lb+1, 0);
#else
  k[0] = int2_as_double(tex3D<int2>(elem.ff_d, indx,   lb  , 0));
  k[1] = int2_as_double(tex3D<int2>(elem.ff_d, indx+1, lb  , 0));
  k[2] = int2_as_double(tex3D<int2>(elem.ff_d, indx,   lb+1, 0));
  k[3] = int2_as_double(tex3D<int2>(elem.ff_d, indx+1, lb+1, 0));
#endif
  
  // Linear interpolation
  //
  double a = (rn - dC*(lb+0)) / dC;
  double b = (dC*(lb+1) - rn) / dC;

  double K = A*(a*k[0] + b*k[2]) + B*(a*k[1] + b*k[3]);

  // Assign the photon energy
  //
  ph = pow(10, K) * hbc;

  // Use the integrated cross section from the differential grid
  //

  xc = 
#if cuREAL == 4
    A*tex1D<float>(elem.ff_0, indx  ) +
    B*tex1D<float>(elem.ff_0, indx+1) ;
#else
    A*int2_as_double(tex1D<int2>(elem.ff_0, indx  )) +
    B*int2_as_double(tex1D<int2>(elem.ff_0, indx+1)) ;
#endif
}


__global__
void testFreeFree
(dArray<cuFP_t> energy,
 dArray<cuFP_t> randsl,
 dArray<cuFP_t> ph, dArray<cuFP_t> xc,
 cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeFreeFree(energy._v[tid], randsl._v[tid], 
		    ph._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}


__device__
void computeColExcite
(cuFP_t E, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem)
{
  if (E < elem.ceEmin or E > elem.ceEmax) {

    xc = 0.0;
    ph= 0.0;

  } else {

    // Interpolate the values
    //
    int indx = std::floor( (E - elem.ceEmin)/elem.ceDelE );
    
    // Sanity check
    //
    if (indx > elem.NColl-2) indx = elem.NColl - 2;
    if (indx < 0)            indx = 0;
    
    double eA   = elem.ceEmin + elem.ceDelE*indx;
    double eB   = elem.ceEmin + elem.ceDelE*(indx+1);
    
    double A = (eB - E)/elem.ceDelE;
    double B = (E - eA)/elem.ceDelE;
    
#if cuREAL == 4
    xc = 
      A*tex3D<float>(elem.ce_d, indx,   0, 0) +
      B*tex3D<float>(elem.ce_d, indx+1, 0, 0) ;
    ph = 
      A*tex3D<float>(elem.ce_d, indx,   1, 0) +
      B*tex3D<float>(elem.ce_d, indx+1, 1, 0) ;
#else
    xc = 
      A*int2_as_double(tex3D<int2>(elem.ce_d, indx  , 0, 0)) +
      B*int2_as_double(tex3D<int2>(elem.ce_d, indx+1, 0, 0)) ;
    ph= 
      A*int2_as_double(tex3D<int2>(elem.ce_d, indx  , 1, 0)) +
      B*int2_as_double(tex3D<int2>(elem.ce_d, indx+1, 1, 0)) ;
#endif
  }
  // DONE
}

__global__ void testColExcite
(dArray<cuFP_t> energy,
 dArray<cuFP_t> ph, dArray<cuFP_t> xc, cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeColExcite(energy._v[tid], ph._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}

__device__
void computeColIonize
(cuFP_t E, cuFP_t& xc, cuIonElement& elem)
{
  if (E < elem.ciEmin or E > elem.ciEmax) {

    xc = 0.0;

  } else {

    // Interpolate the values
    //
    int indx = std::floor( (E - elem.ciEmin)/elem.ciDelE );

    // Sanity check
    //
    if (indx > elem.NIonz-2) indx = elem.NIonz - 2;
    if (indx < 0)            indx = 0;
    
    double eA   = elem.ciEmin + elem.ciDelE*indx;
    double eB   = elem.ciEmin + elem.ciDelE*(indx+1);
    
    double A = (eB - E)/elem.ciDelE;
    double B = (E - eA)/elem.ciDelE;
    
#if cuREAL == 4
    xc = 
      A*tex1D<float>(elem.ci_d, indx  ) +
      B*tex1D<float>(elem.ci_d, indx+1) ;
#else
    xc = 
      A*int2_as_double(tex1D<int2>(elem.ci_d, indx  )) +
      B*int2_as_double(tex1D<int2>(elem.ci_d, indx+1)) ;
#endif
  }
}


__device__
void computePhotoIonize
(cuFP_t rr, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem)
{
  constexpr cuFP_t dC = 1.0/CHCUMK;
  int indx  = rr/dC;
  if (indx > CHCUMK-2) indx = CHCUMK - 2;

  // Linear interpolation
  //
  double a = (rr - dC*(indx+0)) / dC;
  double b = (dC*(indx+1) - rr) / dC;

  ph = a*elem.pi_0[indx+0] + b*elem.pi_0[indx+1];
  xc = elem.piTotl;
}


__global__ void testColIonize
(dArray<cuFP_t> energy, dArray<cuFP_t> xc, cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeColIonize(energy._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}

__device__
void computeRadRecomb
(cuFP_t E, cuFP_t& xc, cuIonElement& elem)
{
  if (E < ionEminGrid or E > ionEmaxGrid) {

    xc = 0.0;

  } else {

    // Interpolate the values
    //
    int indx = std::floor( (E - ionEminGrid)/ionDeltaEGrid );

    // Sanity check
    //
    if (indx > ionRadRecombNumber-2) indx = ionRadRecombNumber - 2;
    if (indx < 0)                    indx = 0;
    
    double eA   = ionEminGrid + ionDeltaEGrid*indx;
    double eB   = ionEminGrid + ionDeltaEGrid*(indx+1);
    
    double A = (eB - E)/ionDeltaEGrid;
    double B = (E - eA)/ionDeltaEGrid;
    
#if cuREAL == 4
    xc = 
      A*tex1D<float>(elem.rc_d, indx  ) +
      B*tex1D<float>(elem.rc_d, indx+1) ;
#else
    xc = 
      A*int2_as_double(tex1D<int2>(elem.rc_d, indx  )) +
      B*int2_as_double(tex1D<int2>(elem.rc_d, indx+1)) ;
#endif
  }
  // DONE
}

__global__
void testRadRecomb
(dArray<cuFP_t> energy, dArray<cuFP_t> xc, cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeRadRecomb(energy._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}


void chdata::testCross(int Nenergy)
{
  // Timers
  //
  Timer serial, cuda;

  // Loop over ions and tabulate statistics
  //
  size_t k = 0;

  thrust::host_vector<cuFP_t> energy_h(Nenergy), randsl_h(Nenergy);

  for (auto v : IonList) {

    IonPtr I = v.second;
    cuIonElement& E = cuIonElem[k];

    // Make an energy grid
    //
    double dE = (I->EmaxGrid - I->EminGrid)/(Nenergy-1) * 0.999;
    for (int i=0; i<Nenergy; i++) {
      energy_h[i] = I->EminGrid + dE*i;
      randsl_h[i] = static_cast<cuFP_t>(rand())/RAND_MAX;
    }

    thrust::device_vector<cuFP_t> energy_d = energy_h;
    thrust::device_vector<cuFP_t> randsl_d = randsl_h;

    // Only free-free for non-neutral species

    thrust::device_vector<cuFP_t> eFF_d(Nenergy), xFF_d(Nenergy);
    thrust::device_vector<cuFP_t> eCE_d(Nenergy), xCE_d(Nenergy);
    thrust::device_vector<cuFP_t> xCI_d(Nenergy), xRC_d(Nenergy);

    unsigned int gridSize  = Nenergy/BLOCK_SIZE;
    if (Nenergy > gridSize*BLOCK_SIZE) gridSize++;

    cuda.start();

    if (E.C>1)
      testFreeFree<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), toKernel(randsl_d),
					     toKernel(eFF_d), toKernel(xFF_d),
					     cuIonElem[k]);

    if (E.C<=E.Z)
      testColExcite<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(eCE_d), toKernel(xCE_d), cuIonElem[k]);
      
    if (E.C<=E.Z)
      testColIonize<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(xCI_d), cuIonElem[k]);
      
    if (E.C>1)
      testRadRecomb<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(xRC_d), cuIonElem[k]);
      
    std::cout << "k=" << k << " delE=" << E.ceDelE << std::endl;

    thrust::host_vector<cuFP_t> eFF_h = eFF_d;
    thrust::host_vector<cuFP_t> xFF_h = xFF_d;
    thrust::host_vector<cuFP_t> eCE_h = eCE_d;
    thrust::host_vector<cuFP_t> xCE_h = xCE_d;
    thrust::host_vector<cuFP_t> xCI_h = xCI_d;
    thrust::host_vector<cuFP_t> xRC_h = xRC_d;
    
    cuda.stop();
    
    std::vector<double> eFF_0(Nenergy, 0), xFF_0(Nenergy, 0);
    std::vector<double> eCE_0(Nenergy, 0), xCE_0(Nenergy, 0);
    std::vector<double> xCI_0(Nenergy, 0), xRC_0(Nenergy, 0);
    
    serial.start();
    
    for (int i=0; i<Nenergy; i++) {
				// Free-free
      auto retFF = I->freeFreeCrossTest(energy_h[i], randsl_h[i], 0);
      if (retFF.first>0.0)
	xFF_0[i]   = (xFF_h[i] - retFF.first )/retFF.first;
      if (retFF.second>0.0)
	eFF_0[i]   = (eFF_h[i] - retFF.second)/retFF.second;

				// Collisional excitation
      auto retCE = I->collExciteCross(energy_h[i], 0).back();
      if (retCE.first>0.0) {
	xCE_0[i]   = (xCE_h[i] - retCE.first )/retCE.first;
	/*
	std::cout << std::setw( 4) << cuZ[k]
		  << std::setw( 4) << cuC[k]
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xCE_h[i]
		  << std::setw(14) << retCE.first
		  << std::endl;
	*/
      }
				// Collisional ionization

      auto retCI = I->directIonCross(energy_h[i], 0);
      if (retCI>0.0) {
	xCI_0[i]   = (xCI_h[i] - retCI)/retCI;
	/*
	std::cout << std::setw( 4) << cuZ[k]
		  << std::setw( 4) << cuC[k]
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xCI_h[i]
		  << std::setw(14) << retCI
		  << std::endl;
	*/
      }

				// Radiative recombination

      auto retRC = I->radRecombCross(energy_h[i], 0).back();
      if (retRC>0.0) {
	xRC_0[i]   = (xRC_h[i] - retRC)/retRC;
	/*
	std::cout << std::setw( 4) << cuZ[k]
		  << std::setw( 4) << cuC[k]
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xRC_h[i]
		  << std::setw(14) << retRC
		  << std::endl;
	*/
      }

      /*
      if (retCE.second>0.0)
	eCE_0[i]   = (eCE_h[i] - retCE.second)/retCE.second;

      if (cuC[k]<=cuZ[k])
	std::cout << std::setw(14) << xCE_h[i]
		  << std::setw(14) << eCE_h[i]
		  << std::endl;
      */
    }

    serial.stop();

    std::sort(xFF_0.begin(), xFF_0.end());
    std::sort(eFF_0.begin(), eFF_0.end());
    std::sort(xCE_0.begin(), xCE_0.end());
    std::sort(eCE_0.begin(), eCE_0.end());
    std::sort(xCI_0.begin(), xCI_0.end());
    std::sort(xRC_0.begin(), xRC_0.end());
    
    std::vector<double> quantiles = {0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99};

    std::cout << "Ion (" << I->Z << ", " << I->C << ")" << std::endl;
    for (auto v : quantiles) {
      int indx = std::min<int>(std::floor(v*Nenergy+0.5), Nenergy-1);
      double FF_xc = 0.0, FF_ph = 0.0, CE_xc = 0.0, CE_ph = 0.0;
      double CI_xc = 0.0, RC_xc = 0.0;
      
      if (E.C>1) {
	FF_xc = xFF_0[indx];
	FF_ph = eFF_0[indx];
	RC_xc = xRC_0[indx];
      }

      if (E.C<=E.Z) {
	CE_xc = xCE_0[indx];
	CE_ph = eCE_0[indx];
	CI_xc = xCI_0[indx];
      }

      std::cout << std::setw(10) << v
		<< " | " << std::setw(14) << FF_xc
		<< " | " << std::setw(14) << FF_ph
		<< " | " << std::setw(14) << CE_xc
		<< " | " << std::setw(14) << CE_ph
		<< " | " << std::setw(14) << CI_xc
		<< " | " << std::setw(14) << RC_xc
		<< std::endl;
    }

    k++;

  } // END: Ion list

  std::cout << std::endl
	    << "Serial time: " << serial() << std::endl
	    << "Cuda time  : " << cuda()   << std::endl;
}


enum cudaInterTypes { 
  nothing    = 0,
  neut_neut  = 1,
  neut_elec  = 2,
  neut_prot  = 3,
  ion_ion    = 4,
  ion_elec   = 5,
  free_free  = 6,
  col_excite = 7,
  col_ionize = 8,
  recombine  = 9,
  elec_elec  = 10
};

const int maxAtomicNumber = 15;
__constant__ cuFP_t cuda_atomic_weights[maxAtomicNumber], cuFloorEV, cuEsu;
__constant__ cuFP_t cuVunit, cuMunit, cuTunit, cuLunit, cuEunit;
__constant__ cuFP_t cuAmu, cuEV, cuLogL, cuCrossfac;
__constant__ bool   cuMeanKE, cuMeanMass, cuNewRecombAlg;

const int coulSelNumT = 2000;
__constant__ cuFP_t coulSelA[coulSelNumT];
__constant__ cuFP_t coulSelTau_i, coulSelTau_f, coulSelTau_z, coulSelDel;

// Coulombic algorithm initialization
//
double cuCA_f(double x)
{
  return 1.0/tanh(x) - 1.0/x;
}

double cuCA_deriv(double x)
{
  if (x>100.0) return 1.0/(x*x);
  double z = sinh(x);
  return 1.0/(x*x) - 1.0/(z*z);
}

double cuCA_func(cuFP_t tau, cuFP_t x)
{
  const double tol = 1.0e-12;
  const int itmax  = 1000;
  double fac       = exp(-2.0*tau), corr;
  
  for (int i=0; i<itmax; i++) {
    corr  = (cuCA_f(x) - fac)/cuCA_deriv(x);
    x    -= corr;
    if (fabs(corr)<tol) break;
  }
  
  return x;
}

// Link static parameters from CollideIon.cc (decide how to link these
// later)
//
static double FloorEv      = 0.05;
static bool   newRecombAlg = false;

void CollideIon::cuda_initialize()
{
  static bool done = false;

  cuda_random_init();

  if (done) return;
  done = true;

  if (c0->cudaDevice>=0) {
    cudaSetDevice(c0->cudaDevice);
  } else {
    std::cerr << "ERROR: c0->cudaDevice not set but CUDA requested"
	      << std::endl;
    MPI_Finalize();
    exit(33);
  }

  thrust::host_vector<cuIonElement> elems;

  int minSp = std::numeric_limits<int>::max();

  for (auto s : SpList) {
    speciesKey k = s.first;
    int Z = k.first;
    int C = k.second;
    // Scan
    bool found = false;
    for (auto & E : ch.cuIonElem) {
      if (E.Z == Z and E.C == C) {
	E.I = s.second;
	minSp = std::min<int>(minSp, s.second);
	elems.push_back(E);
	found = true;
	break;
      }
    }
    if (not found) {
      std::cout << "CollideIon::cuda_initialize: [Z, C] = ["
		<< Z << ", " << C << "] not found" << std::endl;
    }
  }

  for (auto & E : elems) E.I -= minSp;

  cuElems = elems;

  // Coulombic velocity selection
  //
  cuFP_t tau_i = 0.0001, tau_f = 4.0, tau_z = 40.0;
  std::vector<cuFP_t> hA(coulSelNumT);

  cuFP_t del = (log(tau_f) - log(tau_i))/(coulSelNumT-1);
  cuFP_t A   = 0.5/tau_i, T;
  for (int i=0; i<coulSelNumT; i++) {
    T = tau_i*exp(del*i);
    hA[i] = A = cuCA_func(T, A);
  }

  cuda_safe_call(cudaMemcpyToSymbol(coulSelTau_i, &tau_i, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelTau_i");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelTau_f, &tau_f, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelTau_f");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelTau_z, &tau_z, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelTau_z");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelDel, &del, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelDel");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelA, &hA[0], coulSelNumT*sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelA");

  cuda_atomic_weights_init();
}


void CollideIon::cuda_atomic_weights_init()
{
  cudaElasticInit();

  std::vector<cuFP_t> weights(maxAtomicNumber);

  weights[0]  = 0.000548579909; // Mass of electron
  weights[1]  = 1.0079;	       // Hydrogen
  weights[2]  = 4.0026;	       // Helium
  weights[3]  = 6.941;	       // Lithum
  weights[4]  = 9.0122;	       // Beryllium
  weights[5]  = 10.811;	       // Boron
  weights[6]  = 12.011;	       // Carbon
  weights[7]  = 14.007;	       // Nitrogen
  weights[8]  = 15.999;	       // Oxygen
  weights[9]  = 18.998;	       // Florine
  weights[10] = 20.180;	       // Neon
  weights[11] = 22.990;	       // Sodium
  weights[12] = 24.305;	       // Magnesium
  weights[13] = 26.982;	       // Aluminium
  weights[14] = 28.085;	       // Silicon

  cuda_safe_call(cudaMemcpyToSymbol(cuda_atomic_weights, &weights[0], sizeof(cuFP_t)*maxAtomicNumber), 
		 __FILE__, __LINE__, "Error copying cuda_atomic_weights");

  cuFP_t v = TreeDSMC::Vunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuVunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuVunit");
  v = TreeDSMC::Lunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuLunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuLunit");
  v = TreeDSMC::Munit;
  cuda_safe_call(cudaMemcpyToSymbol(cuMunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuMunit");
  v = TreeDSMC::Tunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuTunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuTunit");
  v = TreeDSMC::Eunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuEunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuEunit");
  v = amu;
  cuda_safe_call(cudaMemcpyToSymbol(cuAmu, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuAmu");
  v = eV;
  cuda_safe_call(cudaMemcpyToSymbol(cuEV, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuEV");
  v = logL;
  cuda_safe_call(cudaMemcpyToSymbol(cuLogL, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuLogL");
  v = crossfac;
  cuda_safe_call(cudaMemcpyToSymbol(cuCrossfac, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuCrossfac");
  v = FloorEv;
  cuda_safe_call(cudaMemcpyToSymbol(cuFloorEV, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuFloorEV");
  v = esu;
  cuda_safe_call(cudaMemcpyToSymbol(cuEsu, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuEsu");

  cuda_safe_call(cudaMemcpyToSymbol(cuMeanMass, &MeanMass, sizeof(bool)), 
		 __FILE__, __LINE__, "Error copying cuMeanMass");

  cuda_safe_call(cudaMemcpyToSymbol(cuNewRecombAlg, &newRecombAlg, sizeof(bool)), 
		 __FILE__, __LINE__, "Error copying cuNewRecombAlg");
}  



// CURAND initialization
//
__global__ void initCurand(dArray<curandState> state, unsigned long seed)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < state._s) {
    curand_init(seed, tid, 0, &state._v[tid]);
  }
}

__global__ void cellInitKernel(dArray<cudaParticle> in,    // Particles (all active)
			       dArray<cuFP_t> Ivel2,       // Mean-squared ion velocity (per cell)
			       dArray<cuFP_t> Evel2,       // Mean-squared electron velocity (per cell)
			       dArray<cuFP_t> PiProb,      // Relative electron fraction for BN algorithm (per cell)
			       dArray<cuFP_t> ABrate,      // Plasma rate for BN algorithm (per cell)
			       dArray<cuFP_t> volC,        // Cell's volume
			       dArray<cuFP_t> tauC,        // Cell's time step
			       dArray<cuFP_t> selC,        // Prefactor for per cell pair selection computation
			       dArray<int>    cellI,       // Index to beginning of bodies for this cell
			       dArray<int>    cellN,	   // Number of bodies per cell
			       dArray<cuIonElement> elems, // Species array
			       int epos,		   // Electron position in array
			       unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  const cuFP_t dfac = cuMunit/cuAmu / (cuLunit*cuLunit*cuLunit);

  if (epos+3>=DATTRIB_CUDA) {
    printf("cellInit: epos OAB, epos+3=%d/%d\n", epos+3, DATTRIB_CUDA);
  }

  for (int s=0; s<stride; s++) {

    int c = tid*stride + s;

    if (c < cellI._s) {

      int nbods = cellN._v[c];

      cuFP_t massP = 0.0, numbP = 0.0, massE = 0.0;
      cuFP_t evel2 = 0.0, ivel2 = 0.0, numQ2 = 0.0;
      cuFP_t molPP = 0.0, densQ = 0.0, densE = 0.0;
      cuFP_t meanM = 0.0;

      for (size_t i=0; i<nbods; i++) {
	
	// Sanity check
	//
	if (i + cellI._v[c] >= in._s) {
	  printf("Cell init: Wanted %d/%d\n", i + cellI._v[c], in._s);
	}

	// The particle
	cudaParticle & p = in._v[i + cellI._v[c]];

	// Mass per cell
	massP += p.mass;
	
	// Mass-weighted trace fraction
	// [Iterate through all Trace ionization states]
	cuFP_t ee = 0.0;
	for (int k=0; k<elems._s; k++) {
	  cuIonElement& E = elems._v[k];
	  if (E.I >= 5) printf("Crazy I: %d\n", E.I);
	  cuFP_t ff = p.datr[E.I];
	  cuFP_t ww = ff/cuda_atomic_weights[E.Z];
	  cuFP_t qq = E.C - 1;
	  // Mean number
	  numbP += p.mass * ww;
	  // Electron fraction
	  ee += ww * qq;
	  // For computing mean molecular weight
	  molPP += p.mass * ww;
	  // Charge
	  densE += p.mass * ww * qq;
	  // Charge squared
	  numQ2 += p.mass * ww * qq*qq;
	  if (E.C) densQ += p.mass * ww;
	}
	
	cuFP_t eVel2 = 0.0, iVel2 = 0.0;
	for (int l=0; l<3; l++) {
	  cuFP_t ve  = p.datr[epos+l];
	  if (isnan(ve)) {
	    printf("Weird electron\n");
	  }
	  eVel2 += ve*ve;
	  cuFP_t vi  = p.vel[l];
	  iVel2 += vi*vi;
	}
	
	evel2 += p.mass * ee * eVel2;
	ivel2 += p.mass * iVel2;
	massE += p.mass * ee;
      }
  
      if (numbP>0.0) meanM       = massP/numbP;
      if (massP>0.0) Ivel2._v[c] = ivel2/massP;
      if (massE>0.0) Evel2._v[c] = evel2/massE;
      if (densQ>0.0) numQ2      /= densQ;
      if (molPP>0.0) molPP       = massP/molPP;
      
      cuFP_t ddfac = dfac/volC._v[c];

      densQ *= ddfac;
      densE *= ddfac;

      // Compute per channel Coulombic probabilities
      //
      // Ion probabilities
      //
      cuFP_t muii = meanM/2.0;
      cuFP_t muee = cuda_atomic_weights[0]/2.0;
      cuFP_t muie = cuda_atomic_weights[0] * meanM/(cuda_atomic_weights[0] + meanM);
      
      // Ion-Ion
      PiProb._v[c*4 + 0] =
	densQ +
	densE * pow(2.0*Ivel2._v[c], 1.5) * muii*muii /
	(pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie * numQ2);
      //                                               ^
      //                                               |
      // The density is weighted by q^2 for each species
      
      // Ion-Electron
      PiProb._v[c*4 + 1] =
	densQ * pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie * numQ2 /
	(pow(2.0*Ivel2._v[c], 1.5) * muii*muii) +
	densE ;
      
      // Electron-Ion
      PiProb._v[c*4 + 2] =
	densQ +
	densE * pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie /
	(pow(2.0*Evel2._v[c], 1.5) * muee*muee * numQ2);
      
      // Electron-Electron
      PiProb._v[c*4 + 3] =
	densQ * pow(2.0*Evel2._v[c], 1.5) * muee*muee * numQ2 /
	(pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie) +
	densE;
      
      // Rate coefficients
      ABrate._v[c*4 + 0] = 2.0*M_PI * PiProb._v[c*4 + 0] * cuLogL * pow(numQ2*numQ2, 2.0);
      
      ABrate._v[c*4 + 1] = 2.0*M_PI * PiProb._v[c*4 + 1] * cuLogL * pow(numQ2, 2.0);
      
      ABrate._v[c*4 + 2] = 2.0*M_PI * PiProb._v[c*4 + 2] * cuLogL * pow(numQ2, 2.0);
      
      ABrate._v[c*4 + 3] = 2.0*M_PI * PiProb._v[c*4 + 3] * cuLogL ;

      // Selection scaling based on plasma rate
      cuFP_t dens = massP * ddfac;
      cuFP_t tot_pairs =  dens * nbods / molPP * 
	(1.0/PiProb._v[c*4+0] + 1.0/PiProb._v[c*4+1] + 1.0/PiProb._v[c*4+2]);
      // cuFP_t Prob = dens * cuMunit/cuAmu * sqrt(Ivel2._v[c]) * tauC._v[c];
      cuFP_t Prob = dens * sqrt(Ivel2._v[c]) * tauC._v[c];

      if (isnan(Prob) or isnan(tot_pairs) or tot_pairs<=0.0) {
	printf("tot_pairs=%f dens=%f nbods=%d molPP=%f\n",
	       tot_pairs, dens, nbods, molPP);
      }

      selC._v[c] = Prob/tot_pairs;

    } // END: cell

  } // END: stride

}

__global__ void crossSectionKernel(dArray<cudaParticle> in,       // Particle array
				   dArray<curandState> randS,	  // Cuda random number objects
				   dArray<cuFP_t> cross,	  // Cross section for each interaction
				   dArray<cuFP_t> delph,	  // Inelastic energy change for each interaction
				   dArray<unsigned char>  xspcs,  // Ionization state for each interaction
				   dArray<cudaInterTypes> xtype,  // Interaction type for each interaction
				   dArray<cuFP_t> xctot,	  // Total cross section per cell
				   dArray<int>    i1,		  // Location of first particle in pair
				   dArray<int>    i2,		  // Location of second particle in pair
				   dArray<int>    cc,		  // Cell id each pair in i1 and i2
				   dArray<cuFP_t> Ivel2,
				   dArray<cuFP_t> Evel2,
				   dArray<cuFP_t> PiProb,
				   dArray<cuFP_t> ABrate,
				   dArray<int>    flagI, 
				   dArray<cuFP_t> xsc_H,
				   dArray<cuFP_t> xsc_He,
				   dArray<cuFP_t> xsc_pH,
				   dArray<cuFP_t> xsc_pHe,
				   dArray<cuIonElement> elems,
				   int numxc,                      // Number of possible cross sections
				   int epos,			   // Position of electron velocity in datr
				   unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int Nsp = elems._s;

  if (epos+3>=DATTRIB_CUDA) {
    printf("crossSection: epos OAB, epos+3=%d/%d\n", epos+3, DATTRIB_CUDA);
  }

  for (int s=0; s<stride; s++) {

    int n = tid*stride + s;

    if (n < i1._s) {

      // Sanity check
      //
      if (i1._v[n] >= in._s) {
	printf("cross section: i1 wanted %d/%d\n", i1._v[n], in._s);
      }
      if (i2._v[n] >= in._s) {
	printf("cross section: i2 wanted %d/%d\n", i2._v[n], in._s);
      }

      cudaParticle& p1 = in._v[i1._v[n]];
      cudaParticle& p2 = in._v[i2._v[n]];
      int          cid = cc._v[n];
      
      // Zero the cross-section type
      //
      for (int j=0; j<numxc; j++) xtype._v[n*numxc+j] = nothing;

      cuFP_t Eta1=0.0, Eta2=0.0, Mu1=0.0, Mu2=0.0, Sum1=0.0, Sum2=0.0;

      xctot._v[n] = 0.0;	// Total cross section accumulator

      int J = 0;		// Cross section position counter

      for (int k=0; k<Nsp; k++) {
	cuIonElement& E = elems._v[k];

				// Number fraction of ions
	cuFP_t one = p1.datr[E.I] / cuda_atomic_weights[E.Z];
	cuFP_t two = p2.datr[E.I] / cuda_atomic_weights[E.Z];

				// Electron number fraction
	Eta1 += one * (E.C - 1);
	Eta2 += two * (E.C - 1);

	Sum1 += one;
	Sum2 += two;
      }
				// The number of electrons per particle
      Eta1 /= Sum1;
      Eta2 /= Sum2;
				// The molecular weight
      Mu1 = 1.0/Sum1;
      Mu2 = 1.0/Sum2;


      // Number of atoms in each super particle
      //
      // cuFP_t N1 = p1.mass*cuMunit/(Mu1*cuAmu);
      // cuFP_t N2 = p2.mass*cuMunit/(Mu2*cuAmu);
      // ^
      // |
      // +--- Not used in this implementation; originally for internal energy tracking
      //

      cuFP_t vel   = 0.0;
      cuFP_t eVel0 = 0.0, eVel1 = 0.0, eVel2 = 0.0;
      cuFP_t gVel0 = 0.0, gVel1 = 0.0, gVel2 = 0.0;
      cuFP_t sVel1 = 0.0, sVel2 = 0.0;
      for (int k=0; k<3; k++) {
	cuFP_t del = p1.vel[k] - p2.vel[k];
	vel += del*del;
	del = p1.vel[k] - p2.vel[k];

	cuFP_t rvel0 = p1.datr[epos+k] - p2.datr[epos+k];
	cuFP_t rvel1 = p1.datr[epos+k] - p2.vel[k];
	cuFP_t rvel2 = p2.datr[epos+k] - p1.vel[k];

	eVel0 += rvel0*rvel0;
	eVel1 += rvel1*rvel1;
	eVel2 += rvel2*rvel2;

	rvel1 = p1.datr[epos+k] - p1.vel[k];
	rvel2 = p2.datr[epos+k] - p2.vel[k];

	sVel1 += rvel1*rvel1;
	sVel2 += rvel2*rvel2;

	// Scaled electron relative velocity
	if (cuMeanMass) {
	  rvel0 = p1.datr[epos+k]*sqrt(Eta1/Mu1) - p2.datr[epos+k]*sqrt(Eta2/Mu2);
	  rvel1 = p1.datr[epos+k]*sqrt(Eta1/Mu1) - p2.vel[k];
	  rvel2 = p2.datr[epos+k]*sqrt(Eta2/Mu2) - p1.vel[k];
	  
	  gVel0 += rvel0*rvel0;
	  gVel1 += rvel1*rvel1;
	  gVel2 += rvel2*rvel2;
	}
      }

      // Energy available in the center of mass of the atomic collision
      //
      vel   = sqrt(vel)   * cuVunit;
      eVel0 = sqrt(eVel0) * cuVunit;
      eVel1 = sqrt(eVel1) * cuVunit;
      eVel2 = sqrt(eVel2) * cuVunit;
      sVel1 = sqrt(sVel1) * cuVunit;
      sVel2 = sqrt(sVel2) * cuVunit;

      // Pick scaled relative velocities for mean-mass algorithm
      if (cuMeanMass) {
	gVel0 = sqrt(gVel0) * cuVunit / vel;
	gVel1 = sqrt(gVel1) * cuVunit / vel;
	gVel2 = sqrt(gVel2) * cuVunit / vel;
      }
      // Pick true relative velocity for all other algorithms
      else {
	gVel0 = eVel0 / vel;
	gVel1 = eVel1 / vel;
	gVel2 = eVel2 / vel;
      }

      cuFP_t  m1  = Mu1 * cuAmu;
      cuFP_t  m2  = Mu2 * cuAmu;
      cuFP_t  me  = cuda_atomic_weights[0] * cuAmu;
      cuFP_t mu0  = m1 * m2 / (m1 + m2);
      cuFP_t mu1  = m1 * me / (m1 + me);
      cuFP_t mu2  = m2 * me / (m2 + me);

      cuFP_t kEi  = 0.5 * mu0 * vel * vel;

      cuFP_t facE = 0.5  * cuAmu * cuVunit * cuVunit / cuEV;
      cuFP_t Eion = facE * Ivel2._v[cid] * 0.5*(Mu1 + Mu2);
      cuFP_t Eelc = facE * Evel2._v[cid] * cuda_atomic_weights[0];
      cuFP_t kEe1 = 0.5  * mu1 * eVel2*eVel2;
      cuFP_t kEe2 = 0.5  * mu2 * eVel1*eVel1;
      // cuFP_t kEee = 0.25 * me  * eVel0*eVel0;

      // Loop through all possible ion states
      //
      for (int k1=0; k1<Nsp; k1++) {

	cuIonElement& elem = elems._v[k1];

	int Z = elem.Z;
	int C = elem.C;
	int P = elem.C - 1;
	int I = elem.I;

	cuFP_t fac1 = p1.datr[I] / cuda_atomic_weights[Z] / Sum1;
	cuFP_t fac2 = p2.datr[I] / cuda_atomic_weights[Z] / Sum2;

	// Loop through all possible ion states for second ion
	//
	for (int kk=0; kk<Nsp; kk++) {
	  
	  cuIonElement& elem2 = elems._v[kk];
	  
	  int ZZ = elem2.Z;
	  int CC = elem2.C;
	  int PP = elem2.C - 1;
	  int II = elem2.I;

	  //--------------------------------------------------
	  // Particle 1 interacts with Particle 2
	  //--------------------------------------------------
	  
	  cuFP_t cfac = p1.datr[I] * p2.datr[II];

	  //-------------------------------
	  // *** Both particles neutral
	  //-------------------------------
    
	  if (P==0 and PP==0) {
	    
	    // Geometric cross sections based on
	    // atomic radius
	    cuFP_t crs = (cudaGeometric(Z) + cudaGeometric(ZZ)) * cfac;
	
	    if (crs>0.0) {

	      // Double counting
	      if (Z == ZZ) crs *= 0.5;

	      cross._v[n*numxc+J] = crs*cuCrossfac;
	      xspcs._v[(n*numxc+J)+0] = Z;
	      xspcs._v[(n*numxc+J)+1] = C;
	      xspcs._v[(n*numxc+J)+2] = I;
	      xspcs._v[(n*numxc+J)+3] = ZZ;
	      xspcs._v[(n*numxc+J)+4] = CC;
	      xspcs._v[(n*numxc+J)+5] = II;
	      xtype._v[n*numxc+J] = neut_neut;
	      xctot._v[n]        += crs*cuCrossfac;
	      
	      J++;		// Increment interaction counter
	    }
	  }

	  // --------------------------------------
	  // *** Neutral atom-proton scattering
	  // --------------------------------------

	  cuFP_t crs1 = 0;

	  // Particle 2 is proton
	  //
	  if (ZZ==1 and CC==2) {
	    // Particle 1 is neutral hydrogen
	    if (Z==1 and P==0) crs1 = cudaElasticInterp(log10(kEi), cuPH_Emin,  cuH_H,   xsc_pH );
	    // Particle 1 is neutral helium
	    if (Z==2 and P==0) crs1 = cudaElasticInterp(log10(kEi), cuPHe_Emin, cuPHe_H, xsc_pHe);
	    crs1 *= cuCrossfac * cfac;
	  }
	  
	  // Particle 1 is proton
	  //
	  if (Z==1 and C==2) {
	    // Particle 2 is neutral hydrogen
	    if (ZZ==1 and PP==0) crs1 = cudaElasticInterp(log10(kEi), cuPH_Emin,  cuPH_H,  xsc_pH );
	    // Particle 2 is neutral helium
	    if (ZZ==2 and PP==0) crs1 = cudaElasticInterp(log10(kEi), cuPHe_Emin, cuPHe_H, xsc_pHe);
	    crs1 *= cuCrossfac * cfac;
	  }

	  if (crs1>0.0) {
	    cross._v[n*numxc+J] = crs1;
	    xspcs._v[(n*numxc+J)*6+0] = Z;
	    xspcs._v[(n*numxc+J)*6+1] = C;
	    xspcs._v[(n*numxc+J)*6+2] = I;
	    xspcs._v[(n*numxc+J)*6+3] = ZZ;
	    xspcs._v[(n*numxc+J)*6+4] = CC;
	    xspcs._v[(n*numxc+J)*6+5] = II;
	    xtype._v[n*numxc+J] = neut_prot;
	    xctot._v[n]        += crs1;
	    
	    J++;
	  }

	  // --------------------------------------
	  // *** Ion-ion scattering
	  // --------------------------------------
	  //
	  if (P>0 and PP>0) {
	    cuFP_t kEc  = cuMeanKE ? Eion : kEi;
	    kEc = 2.0*kEc > cuFloorEV ? 2.0*kEc*cuEV : cuFloorEV*cuEV;
	    cuFP_t afac = cuEsu*cuEsu/kEc * 1.0e7;
	    cuFP_t crs  = 2.0 * ABrate._v[cid*4+0] * afac*afac / PiProb._v[cid*4+0] * cfac;

	    if (cfac>0.0) {

	      cross._v[n*numxc+J] = crs;
	      xspcs._v[(n*numxc+J)*6+0] = Z;
	      xspcs._v[(n*numxc+J)*6+1] = C;
	      xspcs._v[(n*numxc+J)*6+2] = I;
	      xspcs._v[(n*numxc+J)*6+3] = ZZ;
	      xspcs._v[(n*numxc+J)*6+4] = CC;
	      xspcs._v[(n*numxc+J)*6+5] = II;
	      xtype._v[n*numxc+J] = ion_ion;
	      xctot._v[n]        += crs;
	      
	      J++;
	    }
	  }

	} // End of inner species loop

	// --------------------------------------
	// *** Neutral atom-electron scattering
	// --------------------------------------
    
	// Particle 1 ION, Particle 2 ELECTRON
	//
	if (Z<=2 and P==0 and Eta2>0.0) {
	  cuFP_t crs = 0.0;

	  // Hydrogen
	  //
	  if (Z==1) crs = cudaElasticInterp(kEe1, cuH_Emin, cuH_H, xsc_H) * gVel2 * Eta2 *
		      cuCrossfac * fac1;
	  // Helium
	  //
	  if (Z==2) crs = cudaElasticInterp(kEe1, cuHe_Emin, cuHe_H, xsc_He) * gVel2 * Eta2 *
		      cuCrossfac * fac1;

	  if (crs>0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspcs._v[(n*numxc+J)*6+0] = Z;
	    xspcs._v[(n*numxc+J)*6+1] = C;
	    xspcs._v[(n*numxc+J)*6+2] = I;
	    xspcs._v[(n*numxc+J)*6+3] = 0;
	    xspcs._v[(n*numxc+J)*6+4] = 0;
	    xspcs._v[(n*numxc+J)*6+5] = 255;
	    xtype._v[n*numxc+J] = neut_elec;
	    xctot._v[n]        += crs;
	    
	    J++;
	  }
	}

	// Particle 2 ION, Particle 1 ELECTRON
	//
	if (P==0 and Eta1>0.0) {
	  cuFP_t crs = 0.0;

	  // Hydrogen
	  //
	  if (Z==1) crs = cudaElasticInterp(kEe2, cuH_Emin, cuH_H, xsc_H) * gVel1 * Eta1 *
		      cuCrossfac * fac1;

	  // Helium
	  //
	  if (Z==2) crs = cudaElasticInterp(kEe2, cuHe_Emin, cuHe_H, xsc_He) * gVel1 * Eta1 *
		      cuCrossfac * fac1;

	  if (crs>0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspcs._v[(n*numxc+J)*6+0] = 0;
	    xspcs._v[(n*numxc+J)*6+1] = 0;
	    xspcs._v[(n*numxc+J)*6+2] = 255;
	    xspcs._v[(n*numxc+J)*6+3] = Z;
	    xspcs._v[(n*numxc+J)*6+4] = C;
	    xspcs._v[(n*numxc+J)*6+5] = I;
	    xtype._v[n*numxc+J] = neut_elec;
	    xctot._v[n]        += crs;
	    
	    J++;
	  }
	}

	// --------------------------------------
	// *** Ion-electron scattering
	// --------------------------------------
	//
	if (P>0 and Eta2>0) {
	  // Particle 1 ION, Particle 2 ELECTRON
	  cuFP_t crs  = 0.0;

	  cuFP_t kEc  = cuMeanKE ? Eelc : kEe1;
	  kEc = 2.0*kEc > cuFloorEV ? 2.0*kEc*cuEV : cuFloorEV*cuEV;
	  cuFP_t afac = cuEsu*cuEsu/kEc * 1.0e7;
	    
	  crs = 2.0 * ABrate._v[cid*4+1] * afac*afac * gVel2 / PiProb._v[cid*4+1];
	  
	  if (fac1*Eta2 > 0.0) {

	    cross._v[n*numxc+J] = crs;
	    xspcs._v[(n*numxc+J)*6+0] = Z;
	    xspcs._v[(n*numxc+J)*6+1] = C;
	    xspcs._v[(n*numxc+J)*6+2] = I;
	    xspcs._v[(n*numxc+J)*6+3] = 0;
	    xspcs._v[(n*numxc+J)*6+4] = 0;
	    xspcs._v[(n*numxc+J)*6+5] = 255;
	    xtype._v[n*numxc+J] = ion_elec;
	    xctot._v[n]        += crs;
	    
	    J++;
	  }
	  
	  // Particle 2 ION, Particle 1 ELECTRON
	  crs  = 0.0;
	  kEc  = cuMeanKE ? Eelc : kEe2;
	  afac = cuEsu*cuEsu/(2.0*kEc > cuFloorEV ? 2.0*kEc : cuFloorEV) * cuEV * 1.0e7;

	  crs = 2.0 * ABrate._v[cid*4+2] * afac*afac * gVel1 / PiProb._v[cid*4+2];

	  if (fac2*Eta1 > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspcs._v[(n*numxc+J)*6+0] = 0;
	    xspcs._v[(n*numxc+J)*6+1] = 0;
	    xspcs._v[(n*numxc+J)*6+2] = 255;
	    xspcs._v[(n*numxc+J)*6+3] = Z;
	    xspcs._v[(n*numxc+J)*6+4] = C;
	    xspcs._v[(n*numxc+J)*6+5] = I;
	    xtype._v[n*numxc+J] = ion_elec;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	} // end: ion-electron scattering


	//-------------------------------
	// *** Free-free
	//-------------------------------
      
	if (Eta1>0.0 and Eta2>0.0) {

	  // Particle 1 ION, Particle 2 ELECTRON
	  cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, ff, ph;
	  cuFP_t rn;
#if cuREAL == 4
	  rn = curand_uniform(&randS._v[n]);
#else
	  rn = curand_uniform_double(&randS._v[n]);
#endif
	  computeFreeFree(ke, rn, ph, ff, elem);

	  cuFP_t crs  = gVel2 * Eta2 * ff * fac1;
	
	  if (crs>0.0) {

	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspcs._v[(n*numxc+J)*6+0] = Z;
	    xspcs._v[(n*numxc+J)*6+1] = C;
	    xspcs._v[(n*numxc+J)*6+2] = I;
	    xspcs._v[(n*numxc+J)*6+3] = 0;
	    xspcs._v[(n*numxc+J)*6+4] = 0;
	    xspcs._v[(n*numxc+J)*6+5] = 255;
	    xtype._v[n*numxc+J] = free_free;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }

	  // Particle 2 ION, Particle 1 ELECTRON
	  ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV;

#if cuREAL == 4
	  rn = curand_uniform(&randS._v[n]);
#else
	  rn = curand_uniform_double(&randS._v[n]);
#endif
	  computeFreeFree(ke, rn, ph, ff, elem);

	  crs = gVel1 * Eta1 * ff * fac2;
	  
	  if (crs>0.0) {

	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspcs._v[(n*numxc+J)*6+0] = 0;
	    xspcs._v[(n*numxc+J)*6+1] = 0;
	    xspcs._v[(n*numxc+J)*6+2] = 255;
	    xspcs._v[(n*numxc+J)*6+3] = Z;
	    xspcs._v[(n*numxc+J)*6+4] = C;
	    xspcs._v[(n*numxc+J)*6+5] = I;
	    xtype._v[n*numxc+J] = free_free;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}
	// end: free-free 

	//-------------------------------
	// *** Collisional excitation
	//-------------------------------
    
	// Particle 1 nucleus has BOUND ELECTRON, Particle 2 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta2>0.0) {
	  cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, ph, xc;
	  computeColExcite(ke, ph, xc, elem);

	  cuFP_t crs = gVel2 * Eta2 * xc * fac1;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspcs._v[(n*numxc+J)*6+0] = Z;
	    xspcs._v[(n*numxc+J)*6+1] = C;
	    xspcs._v[(n*numxc+J)*6+2] = I;
	    xspcs._v[(n*numxc+J)*6+3] = 0;
	    xspcs._v[(n*numxc+J)*6+4] = 0;
	    xspcs._v[(n*numxc+J)*6+5] = 255;
	    xtype._v[n*numxc+J] = col_excite;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}

    
	// Particle 2 nucleus has BOUND ELECTRON, Particle 1 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta1>0) {

	  cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, ph, xc;
	  computeColExcite(ke, ph, xc, elem);

	  cuFP_t crs = gVel1 * Eta1 * xc * fac2;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspcs._v[(n*numxc+J)*6+0] = 0;
	    xspcs._v[(n*numxc+J)*6+1] = 0;
	    xspcs._v[(n*numxc+J)*6+2] = 255;
	    xspcs._v[(n*numxc+J)*6+3] = Z;
	    xspcs._v[(n*numxc+J)*6+4] = C;
	    xspcs._v[(n*numxc+J)*6+5] = I;
	    xtype._v[n*numxc+J] = col_excite;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}	  
	// end: colexcite

	//-------------------------------
	// *** Ionization cross section
	//-------------------------------
	
	// Particle 1 nucleus has BOUND ELECTRON, Particle 2 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta2>0) {
      
	  cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, xc;
	  computeColIonize(ke, xc, elem);

	  cuFP_t crs = gVel2 * Eta2 * xc * fac1;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspcs._v[(n*numxc+J)*6+0] = Z;
	    xspcs._v[(n*numxc+J)*6+1] = C;
	    xspcs._v[(n*numxc+J)*6+2] = I;
	    xspcs._v[(n*numxc+J)*6+3] = 0;
	    xspcs._v[(n*numxc+J)*6+4] = 0;
	    xspcs._v[(n*numxc+J)*6+5] = 255;
	    xtype._v[n*numxc+J] = col_ionize;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}
    
	// Particle 2 nucleus has BOUND ELECTRON, Particle 1 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta1) {
	  
	  cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, xc;
	  computeColIonize(ke, xc, elem);

	  cuFP_t crs = gVel1 * Eta1 * xc * fac2;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspcs._v[(n*numxc+J)*6+0] = 0;
	    xspcs._v[(n*numxc+J)*6+1] = 0;
	    xspcs._v[(n*numxc+J)*6+2] = 255;
	    xspcs._v[(n*numxc+J)*6+3] = Z;
	    xspcs._v[(n*numxc+J)*6+4] = C;
	    xspcs._v[(n*numxc+J)*6+5] = I;
	    xtype._v[n*numxc+J] = col_ionize;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}
	// end: ionize

	//-------------------------------
	// *** Radiative recombination
	//-------------------------------
	
	// The "new" algorithm uses the electron energy of the ion's
	// electron rather than the standard particle partner.
	//
	
	if (cuNewRecombAlg) {
	  
	  // Particle 1 is ION, Particle 2 has ELECTRON
	  //
	  //  +--- Ion charge
	  //  |
	  //  v
	  if (P>0) {
	    cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = Eta1 * xc * fac1;
	
	    if (cuMeanMass) crs *= gVel1;
	    else            crs *= sVel1;

	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspcs._v[(n*numxc+J)*6+0] = Z;
	      xspcs._v[(n*numxc+J)*6+1] = C;
	      xspcs._v[(n*numxc+J)*6+2] = I;
	      xspcs._v[(n*numxc+J)*6+3] = 0;
	      xspcs._v[(n*numxc+J)*6+4] = 0;
	      xspcs._v[(n*numxc+J)*6+5] = 255;
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n]        += crs;
	      
	      J++;
	    }
	  }
      
	  // Particle 2 is ION, Particle 1 has ELECTRON
	  //
	  //  +--- Ion charge
	  //  |
	  //  v
	  if (P>0) {
	    cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = Eta2 * xc * fac2;
	    
	    if (cuMeanMass) crs *= gVel2;
	    else            crs *= sVel2;

	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspcs._v[(n*numxc+J)*6+0] = 0;
	      xspcs._v[(n*numxc+J)*6+1] = 0;
	      xspcs._v[(n*numxc+J)*6+2] = 255;
	      xspcs._v[(n*numxc+J)*6+3] = Z;
	      xspcs._v[(n*numxc+J)*6+4] = C;
	      xspcs._v[(n*numxc+J)*6+5] = I;
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n] += crs;
	      
	      J++;
	    }
	  }
	} // end: new recomb algorithm
	else {
	  // Particle 1 is ION, Particle 2 has ELECTRON
	  //
	  //  +--- Charge of the current subspecies
	  //  |
	  //  |       +--- Electron fraction of partner
	  //  |       |
	  //  V       V
	  if (P>0 and Eta2>0.0) {
	    cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = Eta2 * xc * fac1;

	    if (cuMeanMass) crs *= gVel2;
	    else            crs *= sVel2;

	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspcs._v[(n*numxc+J)*6+0] = Z;
	      xspcs._v[(n*numxc+J)*6+1] = C;
	      xspcs._v[(n*numxc+J)*6+2] = I;
	      xspcs._v[(n*numxc+J)*6+3] = 0;
	      xspcs._v[(n*numxc+J)*6+4] = 0;
	      xspcs._v[(n*numxc+J)*6+5] = 255;
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n]        += crs;
	      
	      J++;
	    }
	  }

	  // Particle 2 is ION, Particle 1 has ELECTRON
	  //
	  //  +--- Charge of the current subspecies
	  //  |
	  //  |       +--- Electron fraction of partner
	  //  |       |
	  //  V       V
	  if (P>0 and Eta1>0.0) {

	    cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = Eta1 * xc * fac2;
	    
	    if (cuMeanMass) crs *= gVel1;
	    else            crs *= sVel1;

	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspcs._v[(n*numxc+J)*6+0] = 0;
	      xspcs._v[(n*numxc+J)*6+1] = 0;
	      xspcs._v[(n*numxc+J)*6+2] = 255;
	      xspcs._v[(n*numxc+J)*6+3] = Z;
	      xspcs._v[(n*numxc+J)*6+4] = C;
	      xspcs._v[(n*numxc+J)*6+5] = I;
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n]        += crs;
	      
	      J++;
	    }
	  }
	} // end: original recomb algorithm
      
      } // end: outer ionization state loop

    } // end: particle loop

  } // end: stride

}

// Return random 3d unit vector
__device__
void cudaUnitVector(cuFP_t *ret, curandState* state)
{
  enum UV {trig, gaus};		// Method choice
  static const UV uv(gaus);	// Pick the method
  
  if (uv == trig) {
#if cuREAL == 4
    cuFP_t cos_th = 1.0 - 2.0*curand_uniform(state);
    cuFP_t phi    = 2.0*M_PI*curand_uniform(state);
#else
    cuFP_t cos_th = 1.0 - 2.0*curand_uniform_double(state);
    cuFP_t phi    = 2.0*M_PI*curand_uniform_double(state);
#endif

    cuFP_t sin_th = sqrt(fabs(1.0 - cos_th*cos_th));

    ret[0] = sin_th*cos(phi);
    ret[1] = sin_th*sin(phi);
    ret[2] = cos_th;
  } else {
    cuFP_t nrm = 0.0;
    for (int i=0; i<3; i++) {
#if cuREAL == 4
      cuFP_t rn = curand_normal(state);
#else
      cuFP_t rn = curand_normal_double(state);
#endif
      ret[i] = rn;
      nrm += ret[i]*ret[i];
    }
    nrm = sqrt(nrm);

    if (nrm>0.0) {
      for (int i=0; i<3; i++) ret[i] /= nrm;
    }
  } 
}



__device__
cuFP_t cuCA_eval(cuFP_t tau)
{
  if      (tau <= coulSelTau_i) return 1.0/(2.0*tau);
  else if (tau >= coulSelTau_z) return 0.0;
  else if (tau >= coulSelTau_f) return 3.0*exp(-2.0*tau);
  else {
    cuFP_t lo  = log(coulSelTau_i);
    cuFP_t hi  = log(coulSelTau_f);
    cuFP_t lv  = log(tau), lo_t, hi_t;
    int indx   = floor( (lv - lo)/coulSelDel );

    if (indx<0) indx = 0;
    if (indx>coulSelNumT-2) indx = coulSelNumT - 2;

    lo_t = log(coulSelTau_i) + coulSelDel*indx;
    hi_t = log(coulSelTau_i) + coulSelDel*(indx+1);

    cuFP_t A = (hi_t - lv)/coulSelDel;
    cuFP_t B = (lv - lo_t)/coulSelDel;
    
    return A*coulSelA[indx-1] + B*coulSelA[indx];
  }
}

__global__ void photoIonizeKernel(dArray<cudaParticle> in,    // Particle array
				  dArray<cuFP_t> dT,	      // Time steps
				  dArray<curandState> randS,  // Cuda random number objects
				  dArray<cuIonElement> elems, // Species map
				  unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int Nsp = elems._s;

  for (int s=0; s<stride; s++) {

    int n = tid*stride + s;

    if (n < in._s) {
      
      // Sanity check
      //
      if (n >= in._s) {
	printf("photoIonize: wanted %d/%d\n", n, in._s);
      }

      cudaParticle& p = in._v[n];

      // Photoionize all subspecies
      //
      for (int s=0; s<Nsp; s++) {
	cuIonElement& elem = elems._v[s];

	int Z = elem.Z;
	int C = elem.C;
	int I = elem.I;
	
	if (C<=Z) {
	  cuFP_t rn, Ep, Pr;
	  // Select random variate and pick a new photon for each body
	  //
#if cuREAL == 4
	  rn = curand_uniform(&randS._v[n]);
#else
	  rn = curand_uniform_double(&randS._v[n]);
#endif
	  computePhotoIonize(rn, Ep, Pr, elem);

	  // Compute the probability and get the residual electron energy
	  Pr *= cuTunit * dT._v[n];
	  cuFP_t ww  = p.datr[I] * Pr;
	  // cuFP_t Gm  = cuMunit/(amu*cuda_atomic_weights[Z]);

	  if (Pr >= 1.0) {	// Limiting case
	    ww = p.datr[I];
	    p.datr[I  ]  = 0.0;
	    p.datr[I+1] += ww;
	  } else {		// Normal case
	    p.datr[I  ] -= ww;
	    p.datr[I+1] += ww;
	  }

	}
      }
    }
  } // End: photoionizing background
}


//! Select tau given random number U in [0,1)
__device__
cuFP_t cuCA_get(cuFP_t tau, cuFP_t U)
{
  cuFP_t A = cuCA_eval(tau);
  if (U<1.0e-14)
    return -1.0;
  else if (A<1.0e-10)
    return 2.0*U - 1.0;
  else if (A>40.0)
    return 1.0 + log(U)/A;
  else
    return log(exp(-A) + 2.0*U*sinh(A))/A;
}

// Return 3d Colombic scattering vector
//
__device__
void cudaCoulombVector(cuFP_t *ret, cuFP_t *rel, cuFP_t W1, cuFP_t W2, cuFP_t Tau,
		       curandState *state)
{
				// Normalize
  cuFP_t rel2 = 0.0;
  for (int i=0; i<3; i++) rel2 += rel[i]*rel[i];
  cuFP_t vfac = sqrt(rel2);
  if (vfac>0.0) for (int i=0; i<3; i++) rel[i] /= vfac;

				// Random generation
#if cuREAL == 4
  cuFP_t rn  = curand_uniform_double(state);
  cuFP_t phi = 2.0*M_PI*curand_uniform(state);
#else
  cuFP_t rn  = curand_uniform_double(state);
  cuFP_t phi = 2.0*M_PI*curand_uniform_double(state);
#endif


  cuFP_t fac    = (W1>W2 ? W1 : W2)/(W1<W2 ? W1 : W2);
  cuFP_t tau    = fac * Tau;
  cuFP_t cosx   = cuCA_get(tau, rn);
  cuFP_t sinx   = sqrt(fabs(1.0 - cosx*cosx));
  cuFP_t cosp   = cos(phi);
  cuFP_t sinp   = sin(phi);
  cuFP_t g_perp = sqrt(rel[1]*rel[1] + rel[2]*rel[2]);

				// Compute randomly-oriented
				// perpendicular vector
  cuFP_t h[3];
  if (g_perp>0.0) {
    h[0] = g_perp * cosp;
    h[1] = -(rel[1]*rel[0]*cosp + rel[2]*sinp)/g_perp;
    h[2] = -(rel[2]*rel[0]*cosp - rel[1]*sinp)/g_perp;
  } else {
    h[0] = 0.0;
    h[1] = cosp;
    h[2] = sinp;
  }
  
  for (int i=0; i<3; i++) rel[i] = rel[i]*cosx - h[i]*sinx;
}

__device__
void cudaDeferredEnergy
(
 cuFP_t E,
 cuFP_t m1,   cuFP_t m2,
 cuFP_t a,    cuFP_t b,
 cuFP_t *E1,  cuFP_t *E2
)
{
  if (m1 < 1.0) {
    E1[1] += a*E/(a + b);
    E2[0] += b*E/(a + b);
  }
  else if (m2 < 1.0) {
    E1[0] += a*E/(a + b);
    E2[1] += b*E/(a + b);
  }
  else {
    E1[0]  += a*E/(a + b);
    E2[0]  += b*E/(a + b);
  }
}

__device__
void cudaScatterTrace
(cuFP_t m1,    cuFP_t m2,
 cuFP_t eta1,  cuFP_t eta2,
 cuFP_t W1,    cuFP_t W2,
 cuFP_t *E1,   cuFP_t *E2,  cuFP_t &totE,
 cuFP_t *v1,   cuFP_t *v2,
 cuFP_t delE,  cuFP_t Tau,
 curandState *state, bool Coulombic
 )
{
  if (cuMeanMass) {

    if (m1<1.0) m1 *= eta1;
    if (m2<1.0) m2 *= eta2;

    // Total effective mass in the collision
    //
    double mt = m1 + m2;

    // Reduced mass (atomic mass units)
    //
    double mu = m1 * m2 / mt;
    
    // Set COM frame
    //
    cuFP_t vcom[3], vrel[3];
    cuFP_t vi = 0.0;
    
    for (size_t k=0; k<3; k++) {
      vcom[k] = (m1*v1[k] + m2*v2[k])/mt;
      vrel[k] = v1[k] - v2[k];
      vi += vrel[k] * vrel[k];
    }

    // Energy in COM
    //
    cuFP_t kE = 0.5*W2*mu*vi;

    // Energy reduced by loss
    //
    cuFP_t vfac = 1.0;
    totE = kE - delE;
    
    // KE is positive
    //
    if (kE>0.0) {
      // More loss energy requested than available?
      //
      if (totE < 0.0) {
	// Add to energy bucket for these particles
	//
	cudaDeferredEnergy(-totE, m1, m2, W1, W2, E1, E2);
	totE = 0.0;
      }
      // Update the outgoing energy in COM
      vfac = sqrt(totE/kE);
      if (isnan(vfac)) {
	printf("totE=%f kE=%f\n", totE, kE);
      }
    }
    // KE is zero (limiting case)
    //
    else {
      if (delE>0.0) {
	// Defer all energy loss
	//
	cudaDeferredEnergy(delE, m1, m2, W1, W2, E1, E2);
	delE = 0.0;
      } else {
	// Apply delE to COM
	//
	vi = -2.0*delE/(W1*mu);
      }
    }
    
    // Assign interaction energy variables
    //
    if (Coulombic)
      cudaCoulombVector(vrel, vrel, W1, W2, Tau, state);
    else
      cudaUnitVector(vrel, state);
    
    vi   = sqrt(vi);
    for (auto & v : vrel) v *= vi;
    //                         ^
    //                         |
    // Velocity in center of mass, computed from v1, v2 and adjusted
    // according to the inelastic energy loss
    //
    bool crazy = false;
    for (size_t k=0; k<3; k++) {
      v1[k] = vcom[k] + m2/mt*vrel[k] * vfac;
      v2[k] = vcom[k] - m1/mt*vrel[k] * vfac;
      if (isnan(v1[k]) or isnan(v2[k])) crazy = true;
    }

    if (crazy) {
      printf("scatter: vcom={%f %f %f} vrel={%f %f %f} m1=%f m2=%f\n",
	     vcom[0], vcom[1], vcom[2],
	     vrel[0], vrel[1], vrel[2], m1, m2);
    }
    
  } // END: MeanMass
  else {

    // Total effective mass in the collision (atomic mass units)
    //
    cuFP_t mt = m1 + m2;

    // Reduced mass (atomic mass units)
    //
    cuFP_t mu = m1 * m2 / mt;

    // Set COM frame
    //
    cuFP_t vcom[3], vrel[3];
    cuFP_t vi = 0.0, vfac = 1.0;
    
    for (size_t k=0; k<3; k++) {
      vcom[k] = (m1*v1[k] + m2*v2[k])/mt;
      vrel[k] = v1[k] - v2[k];
      vi += vrel[k] * vrel[k];
    }
				// Energy in COM
    cuFP_t kE = 0.5*W2*mu*vi;
				// Energy reduced by loss
    cuFP_t totE = kE - delE;

    // KE is positive
    //
    if (kE>0.0) {
      // More loss energy requested than available?
      //
      if (totE < 0.0) {
	cudaDeferredEnergy(-totE, m1, m2, W1, W2, E1, E2);
	totE = 0.0;
      }
      // Update the outgoing energy in COM
      vfac = sqrt(totE/kE);
    }
    // KE is zero (limiting case)
    //
    else {
      if (delE>0.0) {
	// Defer all energy loss
	//
	cudaDeferredEnergy(delE, m1, m2, W1, W2, E1, E2);
	delE = 0.0;
      } else {
	// Apply delE to COM
	//
	vi = -2.0*delE/(W1*mu);
      }
    }

    // Assign interaction energy variables
    //
    if (Coulombic)
      cudaCoulombVector(vrel, vrel, W1, W2, Tau, state);
    else
      cudaUnitVector(vrel, state);
  
    vi   = sqrt(vi);
    for (auto & v : vrel) v *= vi;
    //                         ^
    //                         |
    // Velocity in center of mass, computed from v1, v2 and adjusted
    // according to the inelastic energy loss
    //

    // BEGIN: energy conservation algorithm

    cuFP_t vrat = 1.0, q = W2/W1, cq = 1.0 - W2/W1;

    if (cq > 0.0 and q < 1.0) {

      cuFP_t uu[3];
      cuFP_t v1i2 = 0.0, b1f2 = 0.0, qT = 0.0;
      cuFP_t udif = 0.0, vcm2 = 0.0;

      for (size_t i=0; i<3; i++) {
	uu[i] = vcom[i] + m2/mt*vrel[i]*vfac;
	vcm2 += vcom[i] * vcom[i];
	v1i2 += v1[i] * v1[i];
	b1f2 += uu[i] * uu[i];
	qT   += v1[i] * uu[i];
	udif += (v1[i] - uu[i]) * (v1[i] - uu[i]);
      }
      
      if (v1i2 > 0.0 and b1f2 > 0.0) qT *= q/v1i2;
      
      vrat = 
	( -qT + std::copysign(1.0, qT)*sqrt(qT*qT + cq*(q*b1f2/v1i2 + 1.0)) )/cq;
    }

    // Assign new velocities
    //
    for (int i=0; i<3; i++) {
      double v0 = vcom[i] + m2/mt*vrel[i]*vfac;
    
      v1[i] = cq*v1[i]*vrat + q*v0;
      v2[i] = vcom[i] - m1/mt*vrel[i]*vfac;
    }
  }
    
} // END: cudaScatterTrace

__global__ void partInteractions(dArray<cudaParticle> in,
				 dArray<curandState> randS,
				 dArray<cuFP_t> cross,
				 dArray<cuFP_t> delph,
				 dArray<unsigned char>  xspcs,
				 dArray<cudaInterTypes> xtype,
				 dArray<cuFP_t> xctot,
				 dArray<int>    i1,
				 dArray<int>    i2,
				 dArray<int>    cc,
				 dArray<cuFP_t> selC,
				 dArray<cuFP_t> PiProb,
				 dArray<cuFP_t> ABrate,
				 dArray<int>    flagI, 
				 dArray<cuIonElement> elems,
				 cuFP_t spTau,
				 int numxc,
				 int epos,
				 unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int Nsp = elems._s;

  if (epos+3>=DATTRIB_CUDA) {
    printf("partInteract: epos OAB, epos+3=%d/%d\n", epos+3, DATTRIB_CUDA);
  }

  // Cache new state weights for each particle
  //
  cuFP_t *FF1 = new cuFP_t [Nsp];
  cuFP_t *FF2 = new cuFP_t [Nsp];

  for (int s=0; s<stride; s++) {

    int n = tid*stride + s;

    if (n < i1._s) {

      // Sanity check
      if (i1._v[n] >= in._s) {
	printf("part interactions i1: wanted %d/%d\n", i1._v[n], in._s);
      }
      if (i2._v[n] >= in._s) {
	printf("part interactions i2: wanted %d/%d\n", i2._v[n], in._s);
      }

      /*
      cudaParticle& p1    = in._v[i1._v[n]];
      cudaParticle& p2    = in._v[i2._v[n]];
      */
      cudaParticle  p1    = in._v[i1._v[n]];
      cudaParticle  p2    = in._v[i2._v[n]];
      curandState*  state = &randS._v[n];

      // Cell position in arrays
      //
      int cid = cc._v[n];
      
      cuFP_t PE[3] = {0, 0, 0}, EE[3] = {0, 0, 0};

      // Cross section in system units
      //
      cuFP_t totalXS = xctot._v[n];

      // Upscale rate
      //
      cuFP_t colCf = selC._v[cid] * totalXS;

      // Energy computation
      //
      cuFP_t Eta1=0.0, Eta2=0.0, Mu1=0.0, Mu2=0.0, Sum1=0.0, Sum2=0.0;
      cuFP_t kEe1 = 0.0, kEe2 = 0.0, kEi = 0.0, Mue = cuda_atomic_weights[0];

      // Electron KE
      //
      cuFP_t iE1 = 0.0, iE2 = 0.0;
      for (size_t k=0; k<3; k++) {
	iE1 += p1.datr[epos+k] * p1.datr[epos+k];
	iE2 += p2.datr[epos+k] * p2.datr[epos+k];

	cuFP_t rvel0 = p1.datr[epos+k] - p2.datr[epos+k];
	cuFP_t rvel1 = p1.datr[epos+k] - p2.vel[k];
	cuFP_t rvel2 = p2.datr[epos+k] - p1.vel[k];

	kEi  += rvel0*rvel0;
	kEe1 += rvel1*rvel1;
	kEe2 += rvel2*rvel2;
      }

      // Electron and molecular weight
      //
      for (int k=0; k<Nsp; k++) {
	cuIonElement& E = elems._v[k];

	FF1[k] = p1.datr[E.I];
	FF2[k] = p2.datr[E.I];

				// Number fraction of ions
	cuFP_t one = p1.datr[E.I] / cuda_atomic_weights[E.Z];
	cuFP_t two = p2.datr[E.I] / cuda_atomic_weights[E.Z];

				// Electron number fraction
	Eta1 += one * (E.C - 1);
	Eta2 += two * (E.C - 1);

	Sum1 += one;
	Sum2 += two;
      }
				// The number of electrons per particle
      Eta1 /= Sum1;
      Eta2 /= Sum2;
				// The molecular weight
      Mu1 = 1.0/Sum1;
      Mu2 = 1.0/Sum2;

      cuFP_t m1  = Mu1 * cuAmu;
      cuFP_t m2  = Mu2 * cuAmu;
      cuFP_t me  = Mue * cuAmu;
  
      cuFP_t mu0 = m1 * m2 / (m1 + m2);
      cuFP_t mu1 = m1 * me / (m1 + me);
      cuFP_t mu2 = me * m2 / (me + m2);

				// Finalize the KEs
      iE1 *= 0.5*p1.mass * cuda_atomic_weights[0];
      iE2 *= 0.5*p2.mass * cuda_atomic_weights[0];

      kEi  *= 0.5 * mu0;
      kEe1 *= 0.5 * mu1;
      kEe2 *= 0.5 * mu2;

      cuFP_t maxP = 0.0;
      cudaInterTypes maxT = nothing;

      struct Itype
      {
	unsigned char Z1, Z2, C1, C2, P1, P2, I1, I2;
	cuFP_t W1, W2, N1, N2;
      } IT;

      for (int J=0; J<numxc; J++) {
	cudaInterTypes T = xtype._v[n*numxc+J];
	if (T == nothing) break;
    
	cuFP_t XS          = cross._v[n*numxc+J];
	cuFP_t XE          = delph._v[n*numxc+J];
	double Prob        = XS/totalXS;

	if (Prob < 1.0e-14) continue;

	// Atomic number and array loc for each member of the pair
	//
	unsigned char *u = &(xspcs._v[(n*numxc+J)*6]);

	IT.Z1 = *(u+0);
	IT.C1 = *(u+1);
	IT.I1 = *(u+2);

	IT.Z2 = *(u+3);
	IT.C2 = *(u+4);
	IT.I2 = *(u+5);

	// Traditional ionization state (e.g. C1=1 is neutral)
	//
	IT.P1 = IT.C1 - 1;
	IT.P2 = IT.C2 - 1;

	// Energy loss
	//
	cuFP_t dE = 0.0;

	// Compute weights and number interacting atoms (255 is
	// electron)
	//
	if (IT.I1<255) IT.W1 = p1.datr[IT.I1] / cuda_atomic_weights[IT.Z1];
	else           IT.W1 = Eta1;

	if (IT.I2<255) IT.W2 = p2.datr[IT.I2] / cuda_atomic_weights[IT.Z2];
	else           IT.W2 = Eta2;

	IT.N1 = IT.W1 * cuMunit / cuAmu;
	IT.N2 = IT.W2 * cuMunit / cuAmu;

	// Number of particles in active partition
	//
	cuFP_t N0 = IT.N1 > IT.N2 ? IT.N2 : IT.N1;

	// Select the maximum probability channel
	//
	if (Prob > maxP) {
	  maxT  = T;
	  maxP  = Prob;
	}

	//-----------------------------
	// Parse each interaction type
	//-----------------------------

	if (T == neut_neut) {
	  PE[0] += Prob;
	}

	if (T == neut_elec) {
	  if (IT.I1<255) PE[1] += Prob;
	  else           PE[2] += Prob;
	}
	
	if (T == neut_prot) {
	  PE[0] += Prob;
	}

	if (T == ion_elec) {
	  if (IT.I1<255) PE[1] += Prob;
	  else           PE[2] += Prob;
	}
      
	// Upscale inelastic cross sections by ratio of cross-section
	// and plasma rate scattering
	//
	Prob *= colCf;
	
	if (T == free_free) {
	  
	  dE = XE * Prob;
				// Sanity
	  if (isnan(dE)) {
	    printf("Crazy dE value\n");
	  }

	  if (IT.I1<255) {	// Ion is p1
	    PE[1] += Prob;
	    EE[1] += dE;
	  } else {		// Ion is p2
	    PE[2] += Prob;
	    EE[2] += dE;
	  }
	}

	if (T == col_excite) {
	  
	  dE = XE * Prob;
				// Sanity
	  if (isnan(dE)) {
	    printf("Crazy dE value\n");
	  }

	  if (IT.I1<255) {	// Ion is p1
	    PE[1] += Prob;
	    EE[1] += dE;
	  } else {		// Ion is p2
	    PE[2] += Prob;
	    EE[2] += dE;
	  }
	} // END: col_excite

	if (T == col_ionize) {

	  cuFP_t WW = Prob;	    

	  dE = XE * Prob;
				// Sanity
	  if (isnan(dE)) {
	    printf("Crazy dE value\n");
	  }

	  if (IT.I1<255) {	// Ion is p1

	    if (IT.I1>Nsp-2) {
	      printf("Crazy ionize I1=%d\n", IT.I1);
	    }
	      
	    if (WW < FF1[IT.I1]) {
	      FF1[IT.I1]   -= WW;
	      FF1[IT.I1+1] += WW;
	    } else {
	      WW = FF1[IT.I1];
	      FF1[IT.I1]    = 0.0;
	      FF1[IT.I1+1] += WW;
	    }

	    Prob = WW;
	  
	    // The kinetic energy of the ionized electron is lost
	    // from the COM KE
	    //
	    cuFP_t Echg = iE1 * Prob / cuda_atomic_weights[IT.Z1];

	    // Energy for ionized electron comes from COM
	    dE += Echg * cuEunit / (N0*eV);
	    
				// Sanity
	    if (isnan(dE)) {
	      printf("Crazy dE value\n");
	    }

	    PE[1] += Prob;
	    EE[1] += dE;
	    
	  } // END: ion-electron
	  else {		// Ion is p2

	    if (IT.I2 > Nsp-2) {
	      printf("Crazy ionize I2=%d\n", IT.I2);
	    }

	    if (WW < FF2[IT.I2]) {
	      FF2[IT.I2]   -= WW;
	      FF2[IT.I2+1] += WW;
	    } else {
	      WW = FF2[IT.I2];
	      FF2[IT.I2]    = 0.0;
	      FF2[IT.I2+1] += WW;
	    }

	    Prob = WW;
	  
	    // The kinetic energy of the ionized electron is lost
	    // from the COM KE
	    //
	    cuFP_t Echg = iE2 * Prob / cuda_atomic_weights[IT.Z2];

	    // Energy for ionized electron comes from COM
	    dE += Echg * cuEunit / (N0*eV);

				// Sanity
	    if (isnan(dE)) {
	      printf("Crazy dE value\n");
	    }

	    PE[2] += Prob;
	    EE[2] += dE;
	  } // END: electron-ion
	  
	} // END: ionize
	
	if (T == recombine) {

	  cuFP_t WW = Prob;

	  if (IT.I1<255) {		// Ion is p1

	    if (IT.I1>Nsp-1 or IT.I1==0) {
	      printf("Crazy recombine I1=%d\n", IT.I1);
	    }

	    if (WW < FF1[IT.I1]) {
	      FF1[IT.I1]   -= WW;
	      FF1[IT.I1-1] += WW;
	    } else {
	      WW = FF1[IT.I1];
	      FF1[IT.I1]    = 0.0;
	      FF1[IT.I1-1] += WW;
	    }

	    Prob = WW;		// Update to truncated value

	    // Electron KE lost in recombination is radiated by does not
	    // change COM energy
	    //
	    // cuFP_t Echg = iE1 * Prob / cuda_atomic_weights[IT.Z1];

	    // Electron KE fraction in recombination
	    //
	    // cuFP_t eE = Echg * cuEunit / (N0*cuEV);
	    
				// Sanity
	    if (isnan(dE)) {
	      printf("Crazy dE value\n");
	    }

	    PE[1] += Prob;
	    EE[1] += dE;
	  } // END: ion-electron
	  else {		// Ion is p2
	    cuFP_t WW = Prob;

	    if (IT.I2 > Nsp-1 or IT.I2==0) {
	      printf("Crazy recombine I2=%d\n", IT.I2);
	    }

	    if (WW < FF2[IT.I2]) {
	      FF2[IT.I2]   -= WW;
	      FF2[IT.I2-1] += WW;
	    } else {
	      WW = FF2[IT.I2];
	      FF2[IT.I2]    = 0.0;
	      FF2[IT.I2-1] += WW;
	    }
	  
	    Prob = WW;		// Update to truncated value

	    // Electron KE lost in recombination is radiated by does not
	    // change COM energy
	    //
	    // cuFP_t Echg = iE2 * Prob / cuda_atomic_weights[IT.Z2];
	    
	    // Electron KE radiated in recombination
	    // cuFP_t eE = Echg * cuEunit / (N0*eV);
	    
				// Sanity
	    if (isnan(dE)) {
	      printf("Crazy dE value\n");
	    }

	    PE[2] += Prob;
	    EE[2] += dE;
	  } // END: electro-ion

	} // END: recomb

      } // END: interaction loop

      // Convert energy loss from eV to system units and
      // total energy change for all interation
      //
      cuFP_t totalDE = 0.0;
      for (int i=0; i<3; i++) {
	EE[i] *= cuEV / cuEunit;
	totalDE += EE[i];
      }

      // Normalize probabilities and sum inelastic energy changes
      //
      cuFP_t probTot = 0.0;
      for (int i=0; i<3; i++) probTot += PE[i];

      if (probTot > 0.0) {
	for (int i=0; i<3; i++) PE[i] /= probTot;
      }

      //
      // Select interaction
      //
      cuFP_t Pr;
#if cuREAL == 4
      Pr = curand_uniform(state);
#else
      Pr = curand_uniform_double(state);
#endif
      
      unsigned int J = 2;
      if      (Pr < PE[0])         J = 0;
      else if (Pr < PE[0] + PE[1]) J = 1;

      // Deferred energy
      //
      cuFP_t E1[2] = {0.0, 0.0};
      cuFP_t E2[2] = {0.0, 0.0};
      cuFP_t totE  = 0.0;

      // Parse for Coulombic interaction
      //
      bool Coulombic = false;
      cuFP_t Tau = 0.0;
      if (maxT==ion_elec or maxT==ion_ion) {
    	Coulombic = true;

	cuFP_t dT  = spTau * cuTunit;

	if (maxT == ion_elec) {
	  if (IT.I1 < 255) {
	    cuFP_t pVel = sqrt(2.0 * kEe1 / mu1);
	    cuFP_t afac = esu*esu/(2.0*kEe1 > cuFloorEV*cuEV ? 2.0*kEe1 : cuFloorEV*cuEV);
	    Tau = ABrate._v[cid*4+1]*afac*afac*pVel * dT;
	  } else {
	    cuFP_t pVel = sqrt(2.0 * kEe2 / mu2);
	    cuFP_t afac = esu*esu/(2.0*kEe2 > cuFloorEV*cuEV ? 2.0*kEe2 : cuFloorEV*cuEV);
	    Tau = ABrate._v[cid*4+2]*afac*afac*pVel * dT;
	  }
	} else {
	  cuFP_t pVel = sqrt(2.0 * kEi / mu0);
	  cuFP_t afac = esu*esu/(2.0*kEi > cuFloorEV*cuEV ? 2.0*kEi : cuFloorEV*cuEV);
	  Tau = ABrate._v[cid*4+0]*afac*afac*pVel * dT;
	}
      }

      
      //
      // Apply neutral-neutral scattering and energy loss
      //
      if (J==0) {
	cuFP_t v1[3], v2[3];

	for (int k=0; k<3; k++) {
	  // Both particles are neutrals or ions
	  v1[k]  = p1.vel[k];
	  v2[k]  = p2.vel[k];
	}

	if (IT.W1 >= IT.W2)
	  cudaScatterTrace(Mu1, Mu2, Eta1, Eta2, IT.W1, IT.W2,
			   &E1[0], &E2[0], totE,
			   &v1[0], &v2[0], totalDE, Tau, state, Coulombic);
	else
	  cudaScatterTrace(Mu2, Mu1, Eta2, Eta1, IT.W2, IT.W1,
			   &E2[0], &E1[0], totE,
			   &v2[0], &v1[0], totalDE, Tau, state, Coulombic);

	// Time-step computation
	//
	{
	  cuFP_t dt = totE/totalDE;
	  p1.dtreq = dt < p1.dtreq ? dt : p1.dtreq;
	  p2.dtreq = dt < p2.dtreq ? dt : p2.dtreq;
	}
      
	for (int k=0; k<3; k++) {
	  // Both particles are ions
	  p1.vel[k] = v1[k];
	  p2.vel[k] = v2[k];
				// Sanity
	  if (isnan(v1[k])) {
	    printf("Strange value for p1.vel for J=0\n");
	  }
	  if (isnan(v2[k])) {
	    printf("Strange value for p2.vel for J=0\n");
	  }
	}
      
      } // END: Atom-atom interaction


      // Apply ion/neutral-electron scattering and energy loss
      // Ion is Particle 1, Electron is Particle 2
      //
      if (J==1) {
	cuFP_t v1[3], v2[3];
      
	for (int k=0; k<3; k++) {
	  // Particle 1 is the ion
	  v1[k]  = p1.vel[k];
	  // Particle 2 is the electron
	  v2[k]  = p2.datr[epos+k];
	}

	PE[1] = totalDE;
      
	if (IT.W1 >= IT.W2)
	  cudaScatterTrace(Mu1, Mue*Eta2, Eta1, Eta2, IT.W1, IT.W2,
			   &E1[0], &E2[0], totE,
			   &v1[0], &v2[0], totalDE, Tau, state, Coulombic);
	else
	  cudaScatterTrace(Mue*Eta2, Mu1, Eta2, Eta1, IT.W2, IT.W1,
			   &E2[0], &E1[0], totE,
			   &v2[0], &v1[0], totalDE, Tau, state, Coulombic);
			   

	// Time-step computation
	//
	{
	  cuFP_t dt = totE/totalDE;
	  p1.dtreq = dt < p1.dtreq ? dt : p1.dtreq;
	  p2.dtreq = dt < p2.dtreq ? dt : p2.dtreq;
	}

      
	for (int k=0; k<3; k++) {
	  // Particle 1 is the ion
	  p1.vel[k] = v1[k];
				// Sanity
	  if (isnan(v1[k])) {
	    printf("Strange value for p1.vel for J=1\n");
	  }

	  // Particle 2 is the elctron
	  p2.datr[epos+k] = v2[k];
				// Sanity
	  if (isnan(v2[k])) {
	    printf("Strange value for v2 for J=1\n");
	  }
	}
	
      } // END: PE[1] (Ion-electron interaction)

      //
      // Apply ion/neutral-electron scattering and energy loss
      // Ion is Particle 2, Electron is Particle 1
      //
      if (J==2) {
	cuFP_t v1[3], v2[3];

	for (int k=0; k<3; k++) {
	  // Particle 1 is the elctron
	  v1[k]  = p1.datr[epos+k];
	  // Particle 2 is the ion
	  v2[k]  = p2.vel[k];
	}

	PE[2]  = totalDE;
      
	if (IT.W1 >= IT.W2)
	  cudaScatterTrace(Mue*Eta1, Mu2, Eta1, Eta2, IT.W1, IT.W2,
			   &E1[0], &E2[0], totE,
			   &v1[0], &v2[0], totalDE, Tau, state, Coulombic);
	else
	  cudaScatterTrace(Mu2, Mue*Eta1, Eta2, Eta1, IT.W2, IT.W1,
			   &E2[0], &E1[0], totE,
			   &v2[0], &v1[0], totalDE, Tau, state, Coulombic);

	// Time-step computation
	//
	{
	  cuFP_t dt = totE/totalDE;
	  p1.dtreq = dt < p1.dtreq ? dt : p1.dtreq;
	  p2.dtreq = dt < p2.dtreq ? dt : p2.dtreq;
	}
	
	for (int k=0; k<3; k++) {
	  // Particle 1 is the electron
	  p1.datr[epos+k] = v1[k];

				// Sanity
	  if (isnan(v1[k])) {
	    printf("Strange value for v2 for J=2\n");
	  }

	  // Particle 2 is the ion
	  p2.vel[k] = v2[k];
				// Sanity
	  if (isnan(v2[k])) {
	    printf("Strange value for p2.vel for J=2\n");
	  }
	}
      
      } // END: Electron-Ion interaction

      // Reassign weights and update electron fractions
      //
      Eta1 = Eta2 = 0.0;
      Sum1 = Sum2 = 0.0;
      for (int k=0; k<Nsp; k++) {
	cuIonElement& E = elems._v[k];

	p1.datr[E.I] = FF1[k];
	p2.datr[E.I] = FF2[k];

				// Number fraction of ions
	cuFP_t one = p1.datr[E.I] / cuda_atomic_weights[E.Z];
	cuFP_t two = p2.datr[E.I] / cuda_atomic_weights[E.Z];

				// Electron number fraction
	Eta1 += one * (E.C - 1);
	Eta2 += two * (E.C - 1);

	Sum1 += one;
	Sum2 += two;
      }
				// Recomputation of Eta
      Eta1 /= Sum1;
      Eta2 /= Sum2;
				// Mu should be identical but redone
				// here anyway
      Mu1 = 1.0/Sum1;
      Mu2 = 1.0/Sum2;

      // Electron-electron interactions
      {
	cuFP_t  q = (IT.W1 > IT.W2 ? IT.W2/IT.W1 : IT.W1/IT.W2);
	cuFP_t m1 = cuda_atomic_weights[0];
	cuFP_t m2 = cuda_atomic_weights[0];
	if (cuMeanMass) {
	  m1 *= Eta1;
	  m2 *= Eta2;
	}
	cuFP_t mt = m1 + m2;
	cuFP_t mu = m1 * m2 / mt;

	// Calculate pair's relative speed (pre-collision)
	//
	cuFP_t vcom[3], vrel[3], v1[3], v2[3];
	cuFP_t KEcom = 0.0;
	for (int k=0; k<3; k++) {
	  v1[k]   = p1.datr[epos+k];
	  v2[k]   = p2.datr[epos+k];
	  vcom[k] = (m1*v1[k] + m2*v2[k])/mt;
	  vrel[k] = v1[k] - v2[k];
	  KEcom  += vrel[k] * vrel[k];
	}

	// No point in inelastic collsion for zero velocity . . .
	//
	if (KEcom > 0.0) {

	  // Relative speed
	  //
	  cuFP_t vi = sqrt(KEcom);
	  cuFP_t cr = vi * cuVunit;

	  // COM KE
	  //
	  KEcom *= 0.5 * mu;

	  // Kinetic energy in eV
	  //
	  double kEee = 0.5 * mu * amu * cr * cr / eV;


	  double afac = esu*esu/( (2.0*kEee > cuFloorEV ? 2.0*kEee : cuFloorEV)*cuEV );
	  cuFP_t dT   = spTau * cuTunit;
	  cuFP_t Tau  = ABrate._v[cid*4+3] * afac*afac * cr * dT;

	  if (cuMeanMass)
	    cudaCoulombVector(vrel, vrel, 1.0, 1.0, Tau, state);
	  else
	    cudaCoulombVector(vrel, vrel, IT.W1, IT.W2, Tau, state);

	  for (int k=0; k<3; k++) vrel[k] *= vi;

	  if (cuMeanMass) {
	  
	    // Energy deficit correction
	    //
	    KEcom *= 0.5*mu;
	    cuFP_t delE = p1.datr[epos+3] + p2.datr[epos+3];

	    cuFP_t vfac = 0.0;
	    if (KEcom>delE) {
	      vfac = sqrt(1.0 - delE/KEcom);
	      p1.datr[epos+3] = p2.datr[epos+3] = 0.0;
	    } else {
	      p1.datr[epos+3] = p2.datr[epos+3] = 0.5*(delE - KEcom);
	    }

	    // Assign new electron velocities
	    //
	    for (int k=0; k<3; k++) {
	      p1.datr[epos+k] = vcom[k] + m2/mt*vrel[k] * vfac;
	      p2.datr[epos+k] = vcom[k] - m1/mt*vrel[k] * vfac;
				// Sanity
	      if (isnan(v1[k])) {
		printf("Strange value for v1 for elec-elec\n");
	      }
	      if (isnan(v2[k])) {
		printf("Strange value for v2 for elec-elec\n");
	      }
	    }
	  }
	  // Explicit energy conservation using splitting
	  //
	  else if (q < 1.0) {

	    cuFP_t vrat = 1.0;
	    cuFP_t uu[3], vv[3];
	    for (size_t k=0; k<3; k++) {
	      // New velocities in COM
	      uu[k] = vcom[k] + 0.5*vrel[k];
	      vv[k] = vcom[k] - 0.5*vrel[k];
	    }
	    
	    cuFP_t v1i2 = 0.0, b1f2 = 0.0, b2f2 = 0.0;
	    cuFP_t udif = 0.0, vcm2 = 0.0, v1u1 = 0.0;
	    
	    for (size_t k=0; k<3; k++) {
	      // Difference in Particle 1
	      udif += (v1[k] - uu[k]) * (v1[k] - uu[k]);
	      // COM norm
	      vcm2 += vcom[k] * vcom[k];
	      // Normalizations
	      v1i2 += v1[k]*v1[k];
	      b1f2 += uu[k]*uu[k];
	      b2f2 += vv[k]*vv[k];
	      v1u1 += v1[k]*uu[k];
	    }


	    cuFP_t qT = v1u1 * q;
	    if (v1i2 > 0.0) qT /= v1i2;
	    cuFP_t sgn = qT>=0.0 ? 1.0 : -1.0;
	    vrat = ( -qT + sgn*sqrt(qT*qT + (1.0 - q)*(q*b1f2/v1i2 + 1.0) ) )/(1.0 - q);

	    // New velocities in inertial frame
	    //
	    cuFP_t u1[3], u2[3];
	    for (size_t k=0; k<3; k++) {
	      u1[k] = (1.0 - q)*v1[k]*vrat + q*uu[k];
	      u2[k] = vv[k];
	    }

	    // Assign new electron velocities
	    //
	    for (int k=0; k<3; k++) {
	      p1.datr[epos+k] = u1[k];
	      p2.datr[epos+k] = u2[k];
				// Sanity
	      if (isnan(u1[k])) {
		printf("Strange value for u1 for elec-elec\n");
	      }
	      if (isnan(u2[k])) {
		printf("Strange value for u2 for elec-elec\n");
	      }
	    }
	  }
	  // Explicit momentum conservation
	  //
	  else {
	    bool equal = fabs(q - 1.0) < 1.0e-14;

	    double vfac = 1.0;

	    if (equal) {
	      const cuFP_t tol = 0.95; // eps = 0.05, tol = 1 - eps
	      cuFP_t KE0 = 0.5*IT.W1*m1*m2/mt * vi*vi;
	      cuFP_t dKE = p1.datr[epos+3] + p2.datr[epos+3];

	      // Too much KE to be removed, clamp to tol*KE0
	      // 
	      if (dKE/KE0 > tol) {
		// Therefore, remaining excess is:
		// dKE' = dKE - tol*KE0 = dKE*(1 - tol*KE0/dKE);
		//
		cuFP_t ratk = tol*KE0/dKE;
		p1.datr[epos+3] *= (1.0 - ratk);
		p2.datr[epos+3] *= (1.0 - ratk);

		dKE = tol*KE0;
	      } else {
		p1.datr[epos+3] = p2.datr[epos+3] = 0.0;
	      }

	      vfac = sqrt(1.0 - dKE/KE0);
	    }
	  
	    cuFP_t qKEfac  = 0.5*IT.W1*m1*q*(1.0 - q);
	    cuFP_t deltaKE = 0.0;
	    for (int k=0; k<3; k++) {
	      cuFP_t v0 = vcom[k] + m2/mt*vrel[k]*vfac;
	      deltaKE += (v0 - v1[k])*(v0 - v1[k]) * qKEfac;
	      p1.datr[epos+k] = (1.0 - q)*v1[k] + q*v0;
	      p2.datr[epos+k] = vcom[k] - m1/mt*vrel[k]*vfac;
				// Sanity
	      if (isnan(p1.datr[epos+k])) {
		printf("Strange value for p1 for elec-elec\n");
	      }
	      if (isnan(p2.datr[epos+k])) {
		printf("Strange value for p2 for elec-elec\n");
	      }
	    }

	    // Correct energy for conservation
	    if (!equal) {
	      cuFP_t KE1e = 0.0, KE2e = 0.0;
	      for (int k=0; k<3; k++) {
		KE1e += p1.datr[epos+k] * p1.datr[epos+k];
		KE2e += p2.datr[epos+k] * p2.datr[epos+k];
	      }
	  
	      KE1e *= 0.5 * p1.mass * Eta1/Mu1 * cuda_atomic_weights[0];
	      KE2e *= 0.5 * p2.mass * Eta2/Mu2 * cuda_atomic_weights[0];

	      double wght1 = 0.5;
	      double wght2 = 0.5;

	      wght1 = KE1e/(KE1e + KE2e);
	      wght2 = KE2e/(KE1e + KE2e);
	      
	      p1.datr[epos+3] -= deltaKE * wght1;
	      p2.datr[epos+3] -= deltaKE * wght2;
	    }
	  }

	} // end: momentum conservation

      } // end: electron-electron
      
      in._v[i1._v[n]] = p1;
      in._v[i2._v[n]] = p2;

    } // END: interactions with atoms AND electrons
    
  } // END: stride

  delete [] FF1;
  delete [] FF2;
}

// Allocate one generator per particle (overkill, could be tuned to
// save memory)
//
void CollideIon::cuda_random_init()
{
  if (d_randS.size() < Particles().size()) {
    
    std::cout << "Node " << myid
	      << ": CUDA random: size was " << d_randS.size()
	      << ", new size will be " << Particles().size()
	      << std::endl;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, c0->cudaDevice);

    // Compute the cross sections for each interaction en masse
    //
    int N        = Particles().size(); // Number of particles
    int stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    int gridSize = (N+BLOCK_SIZE*stride-1)/(BLOCK_SIZE*stride);

    d_randS.resize(N);
    initCurand<<<gridSize, BLOCK_SIZE>>>(toKernel(d_randS), seed);
  }
}


void * CollideIon::collide_thread_cuda(void * arg)
{
  cuda_initialize();

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, c0->cudaDevice);

  int id = static_cast<int>(((thrd_pass_arguments*)arg)->id);
  
  thread_timing_beg(id);
  
  if (id==0) {
    std::ostringstream sout;
    sout << "Collide::collide: ENTERING cuda collide_thread, T=" << tnow;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

  // Initialize cell loop diagnostics
  //
  pre_cell_loop(id);

  // Start execution timer
  //
  cellTime[id].start();
  
  // Number of cells to process
  //
  size_t Ncells = cellist[id].size();

  // Structures for cell boundaries and counts
  //
  thrust::host_vector<int>    cellI, cellN, pairs, i1, i2, cc;
  thrust::host_vector<char>   flagI;
  thrust::host_vector<cuFP_t> h_volC(Ncells), h_tauC(Ncells);

  size_t Pcount = 0, Npairs = 0, Count = 0, Pindx = 0, Pindx0 = 0;
  typedef std::pair<long int, int> partIndex;
  std::vector<partIndex> bods;

  // Loop over cells to get count of bodies and cells to process
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {
    
    // The current cell
    //
    pCell *c = cellist[id][j];

    // Skip cell if this time has already been computed
    if (c->time >= tnow) {
      continue;
    }

    h_volC[j] = c->Volume();
    h_tauC[j] = dtime / (1<<c->maxplev);

    auto number = c->bods.size();
    if (number>1) {

      cellI.push_back(Pcount);
      cellN.push_back(number);
      pairs.push_back(Npairs);
      flagI.push_back(1);
      Pcount += number;

      for (size_t b=0; b<number; b++)
	bods.push_back(partIndex(c->bods[b], j));

      Pindx0 = Pindx;
      for (size_t b=0; b<number/2; b++) {
	i1.push_back(Pindx++);
	i2.push_back(Pindx++);
	cc.push_back(Count);
	Npairs++;
      }
      if ((number/2)*2 != number) {
	flagI.back() = 2;
	i1.push_back(Pindx0);
	i2.push_back(Pindx-1);
	cc.push_back(Count);
	flagI.push_back(0);
	Npairs++;
      }
      Count++;
    }
  }

  // Prepare for cudaParticle staging
  //
  if (c0->host_particles.capacity()<Pcount) c0->host_particles.reserve(Pcount);
  c0->host_particles.resize(Pcount);

  // Species map info
  //
  int minSp = std::numeric_limits<int>::max(), maxSp = 0;
  int numNeut = 0, numProt = 0, numIon = 0;
  for (auto s : SpList) {
    speciesKey k = s.first;
    minSp = std::min<int>(minSp, s.second);
    maxSp = std::max<int>(maxSp, s.second);
    if (k.second==1) numNeut++;
    else             numIon++;
    if (k.first==1)  numProt++;
  }
  numProt *= numNeut;
  maxSp++;

  int ePos = use_elec - minSp;	// Electron position in cudaParticle datr

  if (use_elec>=0)
    maxSp = std::max<int>(maxSp, use_elec+3);
  else
    throw GenericError("use_elec must be set to use CUDA Trace implementation",  __FILE__, __LINE__);

  // Copy particles to DEVICE
  //
  thrust::host_vector<cuFP_t> h_tauP(Pcount);
  thrust::host_vector<cuFP_t>::iterator pit = h_tauP.begin();
  thrust::host_vector<cudaParticle>::iterator hit = c0->host_particles.begin();
  
  for (auto b : bods) {
    PartPtr h = Particles()[b.first];
    ParticleHtoD(h, *(hit++), minSp, maxSp);
    *(pit++) = h_tauC[b.second];
  }

#ifdef MASS_DEBUG
  std::cout << "Pcount = " << Pcount << ", Ncount = " << pc << std::endl
	    << "Mass[0] = " << c0->host_particles.front().mass << std::endl
	    << "Mass["  << c0->host_particles.size()-1 << "] = "
	    << c0->host_particles.back().mass << std::endl;
#endif

  thrust::device_vector<cuFP_t>       d_tauP(h_tauP);
  thrust::device_vector<cudaParticle> d_part(c0->host_particles);

#ifdef MASS_DEBUG
  c0->host_particles = d_part;		// Check copy

  std::cout << "1st copy back" << std::endl
	    << "Mass[0] = " << c0->host_particles.front().mass << std::endl
	    << "Mass["  << c0->host_particles.size()-1 << "] = "
	    << c0->host_particles.back().mass << std::endl;
#endif

  // Copy cell boundaries and counts to DEVICE
  //
  thrust::device_vector<int>    d_cellI = cellI;
  thrust::device_vector<int>    d_cellN = cellN;
  thrust::device_vector<int>    d_i1    = i1;
  thrust::device_vector<int>    d_i2    = i2;
  thrust::device_vector<int>    d_cc    = cc;
  thrust::device_vector<int>    d_pairs = pairs;

  // Grid size computation
  //
  int N        = cellI.size();	// Number of cells
  int stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  int gridSize = (N+BLOCK_SIZE*stride-1)/(BLOCK_SIZE*stride);

				// These do not need copying back
  thrust::device_vector<cuFP_t> d_Ivel2(N), d_Evel2(N);
  thrust::device_vector<cuFP_t> d_PiProb(N*4), d_ABrate(N*4);
  thrust::device_vector<cuFP_t> d_volC(h_volC), d_tauC(h_tauC), d_selC(N);


  // Initialize per cell info
  //
  cellInitKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),		// Particle array (input)
     toKernel(d_Ivel2),		// Mean squared ion velocity (output)
     toKernel(d_Evel2),		// Mean squared electron velocity (output)
     toKernel(d_PiProb),	// For BN algorithm (output)
     toKernel(d_ABrate),	// For BN algorithm (output)
     toKernel(d_volC),		// Cell volume (input)
     toKernel(d_tauC),		// Cell time step (input)
     toKernel(d_selC),		// Prefactor for inelastic normalization(output)
     toKernel(d_cellI),		// Cell index (input)
     toKernel(d_cellN),		// Cell body count (output)
     toKernel(cuElems),		// Ionization state info (input)
     ePos,			// Electron location in cudaParticle datr
     stride);			// Stride for processing (input)

#ifdef MASS_DEBUG
  c0->host_particles = d_part;		// Check copy

  std::cout << "2nd copy back" << std::endl
	    << "Mass[0] = " << c0->host_particles.front().mass << std::endl
	    << "Mass["  << c0->host_particles.size()-1 << "] = "
	    << c0->host_particles.back().mass << std::endl;
#endif
  
  // Compute the cross sections for each interaction en masse
  //
  N        = i1.size();	// Number of pairs
  stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  gridSize = (N+BLOCK_SIZE*stride-1)/(BLOCK_SIZE*stride);

				// These do not need copying back
  unsigned int totalXCsize =
    numNeut*numNeut + numNeut*2 + numNeut*numProt*2 + numIon*numIon + 
    numIon*5*2;

  thrust::device_vector<cuFP_t>         d_cross(N*totalXCsize);
  thrust::device_vector<cuFP_t>         d_delph(N*totalXCsize);
  thrust::device_vector<unsigned char>  d_xspcs(N*totalXCsize*6);
  thrust::device_vector<cudaInterTypes> d_xtype(N*totalXCsize);
  thrust::device_vector<int>            d_flagI(flagI);

  crossSectionKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),   toKernel(d_randS),
     toKernel(d_cross),  toKernel(d_delph), toKernel(d_xspcs),
     toKernel(d_xtype),  toKernel(d_cross),
     toKernel(d_i1),     toKernel(d_i2),     toKernel(d_cc),
     toKernel(d_Ivel2),  toKernel(d_Evel2),
     toKernel(d_PiProb), toKernel(d_ABrate), toKernel(d_flagI), 
     toKernel(xsc_H),    toKernel(xsc_He),   toKernel(xsc_pH),  toKernel(xsc_pHe),
     toKernel(cuElems),
     totalXCsize, ePos, stride);

#ifdef MASS_DEBUG
  c0->host_particles = d_part;		// Check copy

  std::cout << "3rd copy back" << std::endl
	    << "Mass[0] = " << c0->host_particles.front().mass << std::endl
	    << "Mass["  << c0->host_particles.size()-1 << "] = "
	    << c0->host_particles.back().mass << std::endl;
#endif

  // Do the interactions
  //
  partInteractions<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),   toKernel(d_randS),
     toKernel(d_cross),  toKernel(d_delph),  toKernel(d_xspcs),
     toKernel(d_xtype),  toKernel(d_cross),
     toKernel(d_i1),     toKernel(d_i2),     toKernel(d_cc),    toKernel(d_selC),
     toKernel(d_PiProb), toKernel(d_ABrate), toKernel(d_flagI), 
     toKernel(cuElems),
     spTau[id], totalXCsize, ePos, stride);

  // Photoionization
  //
  if (use_photoIB) {
    N        = Pcount;	// Number of particles
    stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    gridSize = (N+BLOCK_SIZE*stride-1)/(BLOCK_SIZE*stride);

    photoIonizeKernel<<<gridSize, BLOCK_SIZE>>>
      (toKernel(d_part),  toKernel(d_tauP),  toKernel(d_randS),
       toKernel(cuElems), stride);
  }

  // Finally, copy back particles to host
  // 
  c0->host_particles = d_part;

#ifdef MASS_DEBUG
  std::cout << "4th copy back" << std::endl
	    << "Mass[0] = " << c0->host_particles.front().mass << std::endl
	    << "Mass["  << c0->host_particles.size()-1 << "] = "
	    << c0->host_particles.back().mass << std::endl;
#endif

  // Copy particles to HOST
  //
  hit = c0->host_particles.begin();
#ifdef MASS_DEBUG
  long int ibeg=-1, iend;
#endif
  for (auto p : c0->host_particles) {
    long int curr = p.indx;
    PartPtr h = Particles()[curr];
#ifdef MASS_DEBUG
    if (ibeg<0) ibeg = curr; iend = curr;
    double m0 = h->mass, mh = (*hit).mass;
    if (fabs((mh - m0)/m0) > 1.0e-18) {
      std::cout << "Mass mismatch before: m0=" << m0 << " != mh=" << mh
		<< std::endl;
    }
#endif
    ParticleDtoH(p, h, minSp, maxSp);
#ifdef MASS_DEBUG
    double m1 = h->mass;
    if (fabs((m1 - m0)/m0) > 1.0e-18) {
      std::cout << "Mass mismatch after: m0=" << m0 << " != m1=" << m1
		<< std::endl;
    }
#endif
  }

#ifdef MASS_DEBUG
  std::cout << "After step" << std::endl
	    << "Mass[0] = " << c0->host_particles.front().mass << std::endl
	    << "Mass["  << c0->host_particles.size()-1 << "] = "
	    << c0->host_particles.back().mass << std::endl
	    << "mass[" << ibeg << "] = " << Particles()[ibeg]->mass
	    << "Mass[" << iend << "] = " << Particles()[iend]->mass
	    << std::endl;
#endif

  if (id==0) {
    std::ostringstream sout;
    sout << "Collide::collide: AFTER cell loop, T=" << tnow;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

  cellSoFar[id] = cellTime[id].stop();

  // Diagnostics at end of cell loop
  //
  post_cell_loop(id);

  thread_timing_end(id);
  
  return (NULL);
}