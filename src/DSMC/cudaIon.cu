// -*- C++ -*-

#include <iostream>
#include <iomanip>
#include <numeric>

#include <Ion.H>
#include <Elastic.H>
#include <curand.h>
#include <curand_kernel.h>

#include <Timer.H>
#include <cudaUtil.cuH>

__constant__ cuFP_t cuH_H, cuHe_H, cuPH_H, cuPHe_H;
__constant__ cuFP_t cuH_Emin, cuHe_Emin, cuPH_Emin, cuPHe_Emin;

// Charged particle type
enum cudaElasticType { electron, proton };

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
void cudaElasticInit(thrust::device_vector<cuFP_t> & xsc_H,
		     thrust::device_vector<cuFP_t> & xsc_He,
		     thrust::device_vector<cuFP_t> & xsc_pH,
		     thrust::device_vector<cuFP_t> & xsc_pHe)
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
  eV_H.push_back(13.600541506433);	xs_H.push_back(5.9924);

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

  cuda_safe_call(cudaMemcpyToSymbol(cuHe_Emin, &eV_He[0], sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuHe_Emin");

  cuda_safe_call(cudaMemcpyToSymbol(cuHe_H, &dx, sizeof(cuFP_t)), 
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
cuFP_t cudaElasticInterp(cuFP_t E, dArray<cuFP_t> xsc, int Z,
			 cudaElasticType etype = electron,
			 bool pin = true)
{
  // Bohr cross section (pi*a_0^2) in nm
  //
  const cuFP_t b_cross = 0.00879735542978;

  cuFP_t Emin, H;
  bool logV = false;
  int N = xsc._s;

  if (Z==1) {
    if (etype == electron) {
      H    = cuH_H;
      Emin = cuH_Emin;
    } else {
      H    = cuPH_H;
      Emin = cuPH_Emin;
      E    = log10(E);
      logV = true;
    }
  }
  else if (Z==2)
    if (etype == electron) {
      H    = cuHe_H;
      Emin = cuHe_Emin;
    } else {
      H    = cuPHe_H;
      Emin = cuPHe_Emin;
      E    = log10(E);
      logV = true;
    }
  else {
    return 0.0;
  }

  cuFP_t Emax = Emin + N*H, val = 0.0;

  // Enforce return value to grid boundaries for off-grid ordinates.
  // Otherwise, values will be extrapolated.
  if      (pin and E <= Emin) val = xsc._v[0];
  else if (pin and E >= Emax) val = xsc._v[N-1];
  else {

    int indx = 0;
    if (E >= Emax)      indx = xsc._s - 2;
    else if (E <= Emin) indx = 0;
    else                indx = floor( (E - Emin)/H );
    
    // Sanity
    if (indx<0) indx = 0;
    if (indx>xsc._s - 2) indx = xsc._s - 2;

    cuFP_t a = (E - Emin - H*(indx+0))/H;
    cuFP_t b = (Emin + H*(indx+1) - E)/H;
    
    val = a*xsc._v[indx] + b*xsc._v[indx+1];

    if ((logV and val>3.0) or val>80.0) {
      if (pin)
	printf("E=%e a=%e b=%e val=%e [pinned]\n", E, a, b, val);
      else
	printf("E=%e a=%e b=%e val=%e [extrap]\n", E, a, b, val);
    }
  }

  if (logV) val = exp10(val);
  
  return b_cross * val;
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

__global__
void sanityFreeFree(const cuIonElement elem)
{
  for (int i=0; i<ionEgridNumber; i++) {
#if cuREAL == 4
    cuFP_t t1 = tex1D<float>(elem.ff_0, i);
    cuFP_t t2 = tex3D<float>(elem.ff_d, i, 0,          0);
    cuFP_t t3 = tex3D<float>(elem.ff_d, i, CHCUMK/4,   0);
    cuFP_t t4 = tex3D<float>(elem.ff_d, i, CHCUMK/2,   0);
    cuFP_t t5 = tex3D<float>(elem.ff_d, i, 3*CHCUMK/4, 0);
    cuFP_t t6 = tex3D<float>(elem.ff_d, i, CHCUMK-1,   0);
#else
    cuFP_t t1 = int2_as_double(tex1D<int2>(elem.ff_0, i));
    cuFP_t t2 = int2_as_double(tex3D<int2>(elem.ff_d, i, 0,          0));
    cuFP_t t3 = int2_as_double(tex3D<int2>(elem.ff_d, i, CHCUMK/4,   0));
    cuFP_t t4 = int2_as_double(tex3D<int2>(elem.ff_d, i, CHCUMK/2,   0));
    cuFP_t t5 = int2_as_double(tex3D<int2>(elem.ff_d, i, 3*CHCUMK/4, 0));
    cuFP_t t6 = int2_as_double(tex3D<int2>(elem.ff_d, i, CHCUMK-1,   0));
#endif
    printf("%5d %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
	   i, t1, t2, t3, t4, t5, t6);
  }
}

__global__
void sanityRadRecomb(const cuIonElement elem)
{
  for (int i=0; i<ionEgridNumber; i++) {
#if cuREAL == 4
    cuFP_t t1 = tex1D<float>(elem.rc_d, i);
#else
    cuFP_t t1 = int2_as_double(tex1D<int2>(elem.rc_d, i));
#endif
    printf("%5d %13.6e\n", i, t1);
  }
}

void atomicData::cuda_initialize_textures()
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

    E.IPval = 0.0;
    if (E.C<= E.Z) {
      lQ Q(I->Z, I->C);
      E.IPval = I->getIP(Q);
    }

    // The free-free array
    //
    if (E.C>1) {		// Must NOT BE neutral
      cudaTextureDesc texDesc;

      memset(&texDesc, 0, sizeof(texDesc));
      texDesc.readMode       = cudaReadModeElementType;
      texDesc.filterMode     = cudaFilterModePoint;
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

	// Get index of lower bound
	//
	auto ilb = std::distance(I->freeFreeGrid[i].begin(),
				 std::lower_bound(I->freeFreeGrid[i].begin(),
						  I->freeFreeGrid[i].end(),
						  I->freeFreeGrid[i].back()));
	// End points
	//
	h_buffer1[i                              ] = I->kgrid[0];
	h_buffer1[i + (CHCUMK-1)*I->NfreeFreeGrid] = I->kgrid[ilb];

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
	  // If iterator points to the first element, increment the
	  // upper boundary
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

#ifdef SANITY_CHECK
      std::cout << std::string(60, '-') << std::endl
		<< "Free-free Z=" << E.Z << " C=" << E.C << std::endl
		<< std::string(60, '-') << std::endl;

      sanityFreeFree<<<1, 1>>>(E);

      std::cout << std::string(60, '-') << std::endl;
#endif
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

#ifdef SANITY_CHECK
      std::cout << std::string(60, '-') << std::endl
		<< "Rad-recomb Z=" << E.Z << " C=" << E.C << std::endl
		<< std::string(60, '-') << std::endl;

      sanityRadRecomb<<<1, 1>>>(E);

      std::cout << std::string(60, '-') << std::endl;
#endif
    }

    // The collisional excitation array
    //
    if (E.C <= E.Z and I->NcollideGrid>0) {

      E.ceEmin = I->collideEmin;
      E.ceEmax = I->collideEmax;
      E.ceDelE = I->delCollideE;
      E.NColl  = I->NcollideGrid;

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

    // Collisional ionization texture (1-d)
    //
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
	
	std::cout << "Allocating pi_0[" << k << "]" << std::endl;

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

void atomicData::cuda_initialize_grid_constants()
{
  double Emin, Emax, delE;
  int NE, NR;

  for (auto v : IonList) {
    Emin = v.second->EminGrid;
    Emax = v.second->EmaxGrid;
    delE = v.second->DeltaEGrid;

    NE = NR = 0;
    if (v.first.second>1) {
      NE = v.second->NfreeFreeGrid;
      NR = v.second->NradRecombGrid;
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
(cuFP_t E, cuFP_t rr, cuFP_t& ph, cuFP_t& xc, const cuIonElement elem)
{
  // value of h-bar * c in eV*nm
  //
  constexpr double hbc = 197.327;

  // Check neutrality
  //
  if (elem.C==1) {

    xc = 0.0;
    ph = 0.0;

  } else {

    // Energy grid is in log10
    //
    E = log10(E);

    // Enforce minimum and maximum energies
    //
    if (E<ionEminGrid) E = ionEminGrid;
    if (E>ionEmaxGrid) E = ionEmaxGrid;

    int indx = std::floor( (E - ionEminGrid)/ionDeltaEGrid );
    
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

    // Interpolate the cross section array
    //
    cuFP_t k[4];
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
    ph = exp10(K) * hbc;

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
}

__global__
void testElasticE
(dArray<cuFP_t> energy,
 dArray<cuFP_t> xc,
 dArray<cuFP_t> xsc_H,
 dArray<cuFP_t> xsc_He,
 const cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    if (elem.Z==1 and elem.C==1) {
      xc._v[tid] = cudaElasticInterp(energy._v[tid], xsc_H, 1, cudaElasticType::electron);
    } else if (elem.Z==2 and elem.C==1) {
      xc._v[tid] = cudaElasticInterp(energy._v[tid], xsc_He, 2, cudaElasticType::electron);
    }
    else
      xc._v[tid] = 0.0;
  }

  __syncthreads();
}

__global__
void testFreeFree
(dArray<cuFP_t> energy,
 dArray<cuFP_t> randsl,
 dArray<cuFP_t> ph, dArray<cuFP_t> xc,
 const cuIonElement elem)
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
(cuFP_t E, cuFP_t& ph, cuFP_t& xc, const cuIonElement elem)
{
  if (E < elem.ceEmin or E > elem.ceEmax) {

    xc = 0.0;
    ph = 0.0;

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
 dArray<cuFP_t> ph, dArray<cuFP_t> xc, const cuIonElement elem)
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
(cuFP_t E, cuFP_t& xc, const cuIonElement elem)
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
    if (indx < 0)             indx = 0;
    
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
(cuFP_t rr, cuFP_t& ph, cuFP_t& xc, const cuIonElement elem)
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
(dArray<cuFP_t> energy, dArray<cuFP_t> ph, dArray<cuFP_t> xc, const cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeColIonize(energy._v[tid], xc._v[tid], elem);
    ph._v[tid] = elem.IPval;
  }

  __syncthreads();
}

__device__
void computeRadRecomb
(cuFP_t E, cuFP_t& xc, const cuIonElement elem)
{
  // Energy grid is in log10
  //
  E = log10(E);

  // Grid test
  //
  if (E < ionEminGrid or E > ionEmaxGrid or elem.C==1) {

    xc = 0.0;

  } 
  // Interpolate the values
  //
  else {
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
(dArray<cuFP_t> energy, dArray<cuFP_t> xc, const cuIonElement elem)
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


// This is the production version
//
void atomicData::testCross(int Nenergy,
			   thrust::device_vector<cuIonElement>& cuElems,
			   thrust::device_vector<cuFP_t>& x1,
			   thrust::device_vector<cuFP_t>& x2) {}

// This is the standlone version
//
void atomicData::testCross(int Nenergy)
{
  // Bohr cross section (pi*a_0^2) in nm
  //
  const cuFP_t b_cross = 0.00879735542978;

  // Timers
  //
  Timer serial, cuda;

  // Initial header
  //
  std::string separator(10+(14+3)*8, '-');
  std::cout << separator << std::endl
	    << " Cross-section comparison for " << Nenergy << " samples"
	    << std::endl
	    << " Fractional GPU errors relative to CPU"
	    << std::endl << separator << std::endl;

  // Loop over ions and tabulate statistics
  //
  size_t k = 0;

  thrust::host_vector<cuFP_t> energy_h(Nenergy), randsl_h(Nenergy);

  thrust::device_vector<cuFP_t> xsc_H, xsc_He, xsc_pH, xsc_pHe;

  cudaElasticInit(xsc_H, xsc_He, xsc_pH, xsc_pHe);

  for (auto v : IonList) {
    
    IonPtr I = v.second;
    cuIonElement& E = cuIonElem[k];

    // Make an energy grid
    //
    double dE = (I->EmaxGrid - I->EminGrid)/(Nenergy-1) * 0.999;
    for (int i=0; i<Nenergy; i++) {
      energy_h[i] = exp10(I->EminGrid + dE*i);
      randsl_h[i] = static_cast<cuFP_t>(rand())/RAND_MAX;
    }

    thrust::device_vector<cuFP_t> energy_d = energy_h;
    thrust::device_vector<cuFP_t> randsl_d = randsl_h;

    // Only free-free for non-neutral species
    //
    thrust::device_vector<cuFP_t> eFF_d(Nenergy), xFF_d(Nenergy);
    thrust::device_vector<cuFP_t> eCE_d(Nenergy), xCE_d(Nenergy);
    thrust::device_vector<cuFP_t> eCI_d(Nenergy), xCI_d(Nenergy);
    thrust::device_vector<cuFP_t> xRC_d(Nenergy), xEE_d(Nenergy);

    unsigned int gridSize  = Nenergy/BLOCK_SIZE;
    if (Nenergy > gridSize*BLOCK_SIZE) gridSize++;

    testElasticE<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), toKernel(xEE_d),
					   toKernel(xsc_H), toKernel(xsc_He),
					   cuIonElem[k]);
    if (E.C>1)
      testFreeFree<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), toKernel(randsl_d),
					     toKernel(eFF_d), toKernel(xFF_d),
					     cuIonElem[k]);
    if (E.C<=E.Z)
      testColExcite<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(eCE_d), toKernel(xCE_d),
					      cuIonElem[k]);

    if (E.C<=E.Z)
      testColIonize<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(eCI_d), toKernel(xCI_d),
					      cuIonElem[k]);
      
    if (E.C>1)
      testRadRecomb<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(xRC_d), cuIonElem[k]);
      
    thrust::host_vector<cuFP_t> xEE_h = xEE_d;
    thrust::host_vector<cuFP_t> eFF_h = eFF_d;
    thrust::host_vector<cuFP_t> xFF_h = xFF_d;
    thrust::host_vector<cuFP_t> eCE_h = eCE_d;
    thrust::host_vector<cuFP_t> xCE_h = xCE_d;
    thrust::host_vector<cuFP_t> eCI_h = eCI_d;
    thrust::host_vector<cuFP_t> xCI_h = xCI_d;
    thrust::host_vector<cuFP_t> xRC_h = xRC_d;
    
    cuda.stop();
    
    std::vector<double> xEE_0;
    std::vector<double> xFF_0, eFF_0, eFF_1;
    std::vector<double> eCE_0, xCE_0;
    std::vector<double> eCI_0, xCI_0, xRC_0;

    serial.start();
    
    const bool debug = false;

    Elastic elastic;

    int cnt = 0;

    for (int i=0; i<Nenergy; i++) {
				// Neutral-electron
      auto retEE = 0.0;
      if (E.C==1) retEE = elastic(E.Z, energy_h[i]);
      if (retEE>0.0) {
	xEE_0.push_back((xEE_h[i] - retEE)/retEE);
	if (debug) 
	  std::cout << std::setw(12) << "N-e"
		    << std::setw(16) << energy_h[i]
		    << std::setw(16) << xEE_h[i]/b_cross
		    << std::setw(16) << retEE/b_cross
		    << std::setw(4)  << E.Z
		    << std::setw(4)  << E.C
		    << std::endl;
      }

				// Free-free
      auto retFF = I->freeFreeCrossTest(energy_h[i], randsl_h[i], 0);
      if (retFF.first>0.0 and xFF_h[i]>0.0)
	xFF_0.push_back((xFF_h[i] - retFF.first )/retFF.first);
      if (eFF_h[i]>0.0 and retFF.second>0.0) {
	eFF_0.push_back(eFF_h[i]);
	eFF_1.push_back(retFF.second);
      }

      if (debug and retFF.first>0.0)
	std::cout << std::setw(12) << "Free free"
		  << std::setw( 4) << E.Z
		  << std::setw( 4) << E.C
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xFF_h[i]
		  << std::setw(14) << retFF.first
		  << std::endl;

				// Collisional excitation
      auto retCE = I->collExciteCross(energy_h[i], 0).back();
      if (retCE.first>0.0) {
	xCE_0.push_back((xCE_h[i] - retCE.first )/retCE.first);
	if (debug)
	  std::cout << std::setw(12) << "Excite"
		    << std::setw( 4) << E.Z
		    << std::setw( 4) << E.C
		    << std::setw(14) << energy_h[i]
		    << std::setw(14) << xCE_h[i]
		    << std::setw(14) << retCE.first
		    << std::endl;
      }
				// Collisional ionization
      auto retCI = I->directIonCross(energy_h[i], 0);
      if (retCI>0.0) {
	xCI_0.push_back((xCI_h[i] - retCI)/retCI);
	if (debug)
	  std::cout << std::setw(12) << "Ionize"
		    << std::setw( 4) << E.Z
		    << std::setw( 4) << E.C
		    << std::setw(14) << energy_h[i]
		    << std::setw(14) << xCI_h[i]
		    << std::setw(14) << retCI
		    << std::endl;
      }

				// Radiative recombination
      auto retRC = I->radRecombCross(energy_h[i], 0).back();
      if (retRC>0.0) {
	xRC_0.push_back((xRC_h[i] - retRC)/retRC);
	if (debug)
	  std::cout << std::setw(12) << "Rad recomb"
		    << std::setw( 4) << E.Z
		    << std::setw( 4) << E.C
		    << std::setw(14) << energy_h[i]
		    << std::setw(14) << xRC_h[i]
		    << std::setw(14) << retRC
		    << std::endl;
      }

    }

    serial.stop();

    std::sort(xEE_0.begin(), xEE_0.end());
    std::sort(xFF_0.begin(), xFF_0.end());
    std::sort(eFF_0.begin(), eFF_0.end());
    std::sort(eFF_1.begin(), eFF_1.end());
    std::sort(xCE_0.begin(), xCE_0.end());
    std::sort(eCE_0.begin(), eCE_0.end());
    std::sort(xCI_0.begin(), xCI_0.end());
    std::sort(eCI_0.begin(), eCI_0.end());
    std::sort(xRC_0.begin(), xRC_0.end());
    
    std::vector<double> quantiles = {0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99};

    std::cout << "Ion (" << I->Z << ", " << I->C << ")" << std::endl;

    std::cout << std::setw(10) << "Quantile"
	      << " | " << std::setw(14) << "ne xc"
	      << " | " << std::setw(14) << "ff xc"
	      << " | " << std::setw(14) << "ff ph"
	      << " | " << std::setw(14) << "CE xc"
	      << " | " << std::setw(14) << "CE ph"
	      << " | " << std::setw(14) << "CI_xc"
	      << " | " << std::setw(14) << "CI_ph"
	      << " | " << std::setw(14) << "RC_xc"
	      << std::endl << std::setfill('-')
	      <<          std::setw(10) << '-'
	      << " + " << std::setw(14) << '-'
	      << " + " << std::setw(14) << '-'
	      << " + " << std::setw(14) << '-'
	      << " + " << std::setw(14) << '-'
	      << " + " << std::setw(14) << '-'
	      << " + " << std::setw(14) << '-'
	      << " + " << std::setw(14) << '-'
	      << " + " << std::setw(14) << '-'
	      << std::endl << std::setfill(' ');

    auto Q = [](double v, std::vector<double>& V) {
	       if (V.size()==0) return 0.0;
	       int k = std::min<int>(std::floor(v*V.size()+0.5), V.size()-1);
	       return V[k];
	     };

    for (auto v : quantiles) {
      double FF_xc = 0.0, FF_ph = 0.0, CE_xc = 0.0, CE_ph = 0.0;
      double CI_ph = 0.0, CI_xc = 0.0, RC_xc = 0.0, EE_xc = 0.0;
      
      EE_xc = Q(v, xEE_0);

      if (E.C>1) {
	FF_xc = Q(v, xFF_0);
	FF_ph = (Q(v, eFF_0) - Q(v, eFF_1))/Q(v, eFF_1);
	RC_xc = Q(v, xRC_0);
      }

      if (E.C<=E.Z) {
	CE_xc = Q(v, xCE_0);
	CE_ph = Q(v, eCE_0);
	CI_xc = Q(v, xCI_0);
	CI_ph = Q(v, eCI_0);
      }

      std::cout << std::setw(10) << v
		<< " | " << std::setw(14) << EE_xc
		<< " | " << std::setw(14) << FF_xc
		<< " | " << std::setw(14) << FF_ph
		<< " | " << std::setw(14) << CE_xc
		<< " | " << std::setw(14) << CE_ph
		<< " | " << std::setw(14) << CI_xc
		<< " | " << std::setw(14) << CI_ph
		<< " | " << std::setw(14) << RC_xc
		<< std::endl;
    }

    k++;

  } // END: Ion list

  std::cout << separator << std::endl
	    << "Serial time: " << serial() << std::endl
	    << "Cuda time  : " << cuda()   << std::endl
	    << separator << std::endl;
}


void atomicData::testFreeFreeGen(unsigned short ZZ, unsigned short CC,
				 double E, int num,
				 std::vector<double>& xc,
				 std::vector<double>& ph)
{
  int Z = static_cast<int>(ZZ), C = static_cast<int>(CC);

  thrust::host_vector<cuFP_t> energy_h(num), randsl_h(num);

  for (int i=0; i<num; i++) {
    energy_h[i] = E;
    randsl_h[i] = static_cast<cuFP_t>(rand())/RAND_MAX;
  }

  thrust::device_vector<cuFP_t> xsc_H, xsc_He, xsc_pH, xsc_pHe;

  // Loop over ions
  //
  size_t k = 0;

  for (auto v : IonList) {
    
    IonPtr I = v.second;
    cuIonElement& E = cuIonElem[k];

    if (E.Z == Z and E.C == C) {

      thrust::device_vector<cuFP_t> energy_d = energy_h;
      thrust::device_vector<cuFP_t> randsl_d = randsl_h;
      thrust::device_vector<cuFP_t> eFF_d(num), xFF_d(num);

      unsigned int gridSize  = num/BLOCK_SIZE;
      if (num > gridSize*BLOCK_SIZE) gridSize++;
      
      if (E.C>1) {
	testFreeFree<<<gridSize, BLOCK_SIZE>>>
	  (toKernel(energy_d), toKernel(randsl_d),
	   toKernel(eFF_d), toKernel(xFF_d),
	   cuIonElem[k]);
	
	thrust::host_vector<cuFP_t> xFF_h = xFF_d;
	thrust::host_vector<cuFP_t> eFF_h = eFF_d;

	xc.resize(num);
	ph.resize(num);
	thrust::copy(xFF_h.begin(), xFF_h.end(), xc.begin());
	thrust::copy(eFF_h.begin(), eFF_h.end(), ph.begin());
	
	return;
      }
    }
      
    k++;
  } // END: Ion list

  return;
}

void atomicData::testCrossCompare(int Nenergy, double Emin, double Emax,
				  bool eVout, std::string scaling)
{
  // Convert from Rydberg to eV
  //
  const double ryd      = 27.2113845/2.0;
  const double eVerg    = 1.60217733e-12;

  // Bohr cross section (pi*a_0^2) in nm
  //
  const cuFP_t b_cross  = 0.00879735542978;

  // Timers
  //
  Timer serial, cuda;

  // Set energy grid
  //
  Emin = log(Emin);
  Emax = log(Emax);

  double dE = (Emax - Emin)/(Nenergy - 1);

  // Initial header
  //
  int nf = 6;			// Number of output fields

  std::string separator(16*nf, '-'); separator[0] = '#';

  std::cout << separator << std::endl
	    << std::left
	    << std::setw(16) << "# Energy"
	    << std::setw(16) << "Elastic"
	    << std::setw(16) << "Free-free"
	    << std::setw(16) << "Collisional"
	    << std::setw(16) << "Ionization"
	    << std::setw(16) << "Recombination"
	    << std::endl;

  std::cout << std::setw(16) << "# [1]";
  for (int i=1; i<nf; i++) {
    std::ostringstream sout; sout << " [" << i+1  << "]";
    std::cout << std::setw(16) << sout.str();
  }
  std::cout << std::endl << separator << std::endl;


  // Loop over ions and tabulate statistics
  //

  thrust::host_vector<cuFP_t> energy_h(Nenergy), randsl_h(Nenergy);

  thrust::device_vector<cuFP_t> xsc_H, xsc_He, xsc_pH, xsc_pHe;

  cudaElasticInit(xsc_H, xsc_He, xsc_pH, xsc_pHe);

  size_t k = 0;
  for (auto v : IonList) {
    
    unsigned short Z = v.first.first;
    unsigned short C = v.first.second;

    IonPtr I = v.second;
    cuIonElement& E = cuIonElem[k];

    // Ion type
    std::cout << separator << std::endl
	      << "# Q=(" << Z << ", " << C << ")" << std::endl
	      << separator << std::endl;


    // Make an energy grid
    //
    for (int i=0; i<Nenergy; i++) {
      energy_h[i] = exp(Emin + dE*i) * ryd;
      randsl_h[i] = static_cast<cuFP_t>(rand())/RAND_MAX;
    }

    thrust::device_vector<cuFP_t> energy_d = energy_h;
    thrust::device_vector<cuFP_t> randsl_d = randsl_h;

    // Only free-free for non-neutral species
    //
    thrust::device_vector<cuFP_t> eFF_d(Nenergy), xFF_d(Nenergy);
    thrust::device_vector<cuFP_t> eCE_d(Nenergy), xCE_d(Nenergy);
    thrust::device_vector<cuFP_t> eCI_d(Nenergy), xCI_d(Nenergy);
    thrust::device_vector<cuFP_t> xRC_d(Nenergy), xEE_d(Nenergy);

    unsigned int gridSize  = Nenergy/BLOCK_SIZE;
    if (Nenergy > gridSize*BLOCK_SIZE) gridSize++;

    cuda.start();

    testElasticE<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), toKernel(xEE_d),
					   toKernel(xsc_H), toKernel(xsc_He),
					   cuIonElem[k]);
    if (C>1)
      testFreeFree<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), toKernel(randsl_d),
					     toKernel(eFF_d), toKernel(xFF_d),
					     cuIonElem[k]);
    if (C<=Z)
      testColExcite<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(eCE_d), toKernel(xCE_d),
					      cuIonElem[k]);

    if (C<=Z)
      testColIonize<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(eCI_d), toKernel(xCI_d),
					      cuIonElem[k]);
      
    if (C>1)
      testRadRecomb<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(xRC_d), cuIonElem[k]);
      
    thrust::host_vector<cuFP_t> xEE_h = xEE_d;
    thrust::host_vector<cuFP_t> eFF_h = eFF_d;
    thrust::host_vector<cuFP_t> xFF_h = xFF_d;
    thrust::host_vector<cuFP_t> eCE_h = eCE_d;
    thrust::host_vector<cuFP_t> xCE_h = xCE_d;
    thrust::host_vector<cuFP_t> eCI_h = eCI_d;
    thrust::host_vector<cuFP_t> xCI_h = xCI_d;
    thrust::host_vector<cuFP_t> xRC_h = xRC_d;
    
    cuda.stop();
    
    serial.start();
    
    const bool debug = false;

				// No scaling
    double scale = 1.0;
				// Scale to Born xc
    if (scaling.find("born") ==0) scale = 1.0/b_cross;
				// Scale to MB
    if (scaling.find("mbarn")==0) scale = 1.0e+04;

    Elastic elastic;

    for (int i=0; i<Nenergy; i++) {

      double EeV  = energy_h[i];
      double ERyd = EeV / ryd;

      std::cout << std::setw(16) << (eVout ? EeV : ERyd)
		<< std::setw(16) << xEE_h[i]   * scale
		<< std::setw(16) << xFF_h[i]   * scale
		<< std::setw(16) << xCE_h[i]   * scale
		<< std::setw(16) << xCI_h[i]   * scale
		<< std::setw(16) << xRC_h[i]   * scale
		<< std::endl;
    }

    k++;
  } // END: Ion list

  std::cout << separator << std::endl
	    << "Serial time: " << serial() << std::endl
	    << "Cuda time  : " << cuda()   << std::endl
	    << separator << std::endl;
}
