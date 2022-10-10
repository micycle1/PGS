import processing.javafx.*;
import micycle.pgs.*;
import micycle.pgs.PGS_Coloring.ColoringAlgorithm;
import java.util.List;
import org.tinfour.standard.IncrementalTin;

String leafSVG = "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"305\" height=\"343\"><path fill=\"#c00\" d=\"M 139.6101,342.5239 C 138.09029,341.76698 136.09966,342.05783 137.59399,340.22523 C 138.86595,338.83737 142.23737,324.89724 144.24442,313.60056 C 145.00331,309.32928 145.13365,304.62955 145.40645,303.15669 C 146.09197,299.05508 146.47771,294.91742 147.30256,290.83828 C 147.68667,289.07086 148.16486,284.58329 148.36521,280.86591 L 148.72949,274.10701 L 147.14116,273.5536 C 145.97612,273.14768 144.64183,273.20391 142.13456,273.76459 C 134.04268,275.57409 129.40387,276.09142 122.37263,275.96846 C 115.78409,275.85326 114.93974,275.94821 112.91099,277.03251 C 110.5298,278.30515 106.44019,279.83822 105.39222,279.85105 C 105.04203,279.85533 103.96346,279.13581 102.99544,278.25209 C 100.42139,275.90227 98.346014,276.18563 94.895203,279.35804 C 94.37134,279.83965 93.017128,280.51088 91.885826,280.84964 C 89.744431,281.49088 79.151118,281.70765 77.67237,281.1405 C 76.803259,280.80717 76.458281,278.25209 77.282376,278.25209 C 78.432703,278.25209 80.575757,275.53465 80.575757,274.07601 C 80.575757,271.42051 79.012493,270.5466 74.997639,270.95773 C 73.09867,271.15218 70.640904,271.68881 69.536023,272.15023 C 67.182538,273.13308 64.253161,272.9218 60.661636,271.51017 C 58.194498,270.54048 58.116577,270.54312 54.499271,271.71992 C 48.539523,273.65878 41.994029,273.09425 41.994029,270.64137 C 41.994029,270.08489 42.354016,268.93377 42.794012,268.08335 C 44.161514,265.4403 43.334692,264.72221 37.554987,263.53314 C 35.634204,263.13798 34.665429,261.18433 35.5609,259.51196 C 36.782223,257.23113 35.410434,256.72123 27.631286,256.56446 C 22.844518,256.468 20.187364,256.18092 19.138689,255.64693 C 15.583897,253.74395 11.951085,252.01062 8.2350117,250.44643 C 0.88390895,247.42919 -1.2436069,245.40664 0.65842003,243.24365 C 1.4511186,242.34213 2.0715879,242.28028 7.0210827,242.60944 C 11.844923,242.93024 12.922216,242.83341 15.967385,241.80535 C 18.977242,240.78917 19.83622,240.70815 22.544445,241.18495 C 26.777174,241.93015 28.343646,241.5196 30.727147,239.04041 L 32.758366,236.9276 L 36.974277,237.27856 C 40.927621,238.17722 44.714635,235.97593 48.590055,236.61463 C 51.929251,237.77808 56.732946,235.52161 58.403639,232.00485 C 59.122434,230.49191 59.115099,230.10247 58.332903,228.23845 C 57.304753,225.78822 54.17101,222.66925 52.254244,222.18842 C 51.509669,222.00163 50.041533,222.07795 48.99173,222.358 C 47.941986,222.63805 46.509323,222.73181 45.808095,222.56636 C 45.106817,222.40089 42.514984,220.22168 40.048385,217.72365 C 35.614319,213.233 35.563738,213.15118 35.563738,210.46868 C 35.563738,208.20323 35.320562,207.52603 34.090126,206.3644 C 32.623538,204.97986 29.644768,203.27977 26.877978,202.24817 C 24.597522,200.84318 23.859327,198.20804 21.664563,196.66571 C 18.990279,195.08057 17.403553,192.2996 15.225648,190.16386 C 13.856107,189.31787 14.777408,187.75109 17.729708,185.90535 C 19.41514,184.85164 21.065747,183.98952 21.397656,183.98952 C 22.214246,183.98952 23.563752,182.3913 24.026091,180.87675 C 24.279759,180.04562 23.974431,178.72497 23.097446,176.85988 C 22.378811,175.33163 21.030104,172.21539 20.10029,169.93491 C 18.578345,166.20205 18.464101,165.54616 18.954899,163.35733 C 19.595243,160.50157 19.275655,158.83676 17.495326,155.75391 C 16.650437,154.2909 16.314942,153.10147 16.543808,152.38064 C 17.161969,150.43396 14.847605,144.03874 12.775437,141.96763 C 11.079294,140.27235 10.203169,137.4927 11.316194,137.33794 C 11.537237,137.30721 12.635693,137.22503 13.757132,137.15533 C 17.831581,136.90213 18.634722,133.7563 16.05435,128.15755 L 14.636205,125.08048 L 15.722491,122.69721 C 17.043869,119.79809 17.074895,118.70702 15.905702,116.25647 C 15.058895,114.48161 15.058895,114.24591 15.905702,112.47105 C 17.102704,109.96217 17.047456,107.78312 15.757165,106.61606 C 12.143594,101.37877 10.677643,94.974464 7.1671618,89.7318 C 5.2153027,87.26865 5.1531399,85.269737 7.0201734,85.004554 C 8.1808428,84.839701 13.591928,87.483648 13.593057,88.216155 C 16.026437,90.943313 18.950557,92.040853 22.53887,92.119385 C 26.111319,92.135779 26.523187,92.259777 28.165371,93.813455 C 30.388964,95.917234 31.881541,96.689529 33.723714,96.689529 C 34.684025,96.689529 35.768272,97.324323 36.953753,98.580655 C 39.725404,101.03183 43.228517,102.23528 46.207472,104.44132 C 49.034707,106.66139 49.547848,106.8656 52.298831,106.8656 C 53.950987,106.8656 55.681324,107.01076 56.143873,107.18819 C 56.606482,107.36561 57.732688,108.99244 58.646645,110.80337 C 60.284281,114.21842 63.127007,116.74624 64.996667,120.00997 C 66.317824,122.50804 69.220182,125.07541 70.723002,125.07541 C 72.172173,125.07541 75.200565,121.77037 75.797942,119.53677 C 76.541068,116.75845 77.173538,115.97051 78.66076,115.97051 C 79.812586,115.97051 83.255034,119.0083 83.255034,120.02475 C 83.255034,120.31182 83.996022,121.44515 84.901625,122.54326 C 86.298815,124.23742 86.870841,124.53983 88.678239,124.53983 C 90.986548,124.53983 92.034843,125.38047 95.195974,129.76649 C 96.333971,131.34543 97.434626,132.20747 98.772231,132.56747 C 102.87687,133.84224 105.35941,140.28284 110.02679,138.47629 C 111.97763,137.43275 112.25274,135.06981 111.1089,129.1816 C 110.52712,126.18678 110.05042,122.47603 110.04951,120.93549 C 110.04866,119.39496 109.66096,117.16798 109.18807,115.98665 C 108.3637,113.92742 108.37715,113.75586 109.51451,111.8292 C 111.09905,109.14494 110.67781,105.52588 108.20933,100.61602 C 107.19763,98.603734 106.20876,95.833919 106.01177,94.460872 C 105.6879,92.203187 105.77407,91.900986 106.91303,91.301096 C 108.50024,90.465088 108.95765,89.779863 108.968,88.222694 C 108.97603,87.009985 106.48279,83.371219 104.47804,81.669845 C 103.65449,80.970919 103.3835,79.704486 103.04167,74.957695 C 102.80997,71.739823 102.24196,67.789177 101.77936,66.178488 C 100.3633,61.247704 100.88609,57.591993 103.00722,57.591993 C 104.13359,57.591993 106.02458,59.509505 106.52368,61.157769 C 106.71616,61.79339 107.94413,63.299715 109.25259,64.505162 C 111.37158,66.457403 111.92249,66.696893 114.29441,66.696893 C 116.10577,66.696893 117.76307,67.125087 119.47685,68.035851 C 121.85203,69.1231 124.2044,69.524812 126.39156,70.981556 C 127.50415,71.865268 128.94674,72.588304 129.59731,72.588304 C 131.01799,72.588304 132.55391,71.029132 132.55391,69.586989 C 132.55391,69.013653 133.27732,65.58361 134.16148,61.964672 C 135.04565,58.345745 135.76905,54.668611 135.76905,53.793284 C 135.76905,51.945873 137.09739,48.487082 137.80687,48.487082 C 138.89412,48.487082 141.12763,45.790571 141.12763,44.47798 C 141.12763,43.717697 141.73047,41.90661 142.46727,40.453338 C 143.20407,39.000076 143.80691,37.216185 143.80691,36.489139 C 143.80691,34.414907 144.7318,32.493271 146.51251,30.867668 C 147.83327,29.661961 148.24089,28.741412 148.61012,26.130572 C 149.10135,23.171811 150.09379,20.408645 150.46183,17.423285 C 150.73902,14.783822 151.2208,13.387188 152.45413,11.64765 C 153.34882,10.385745 154.48313,8.216643 154.97482,6.8273997 C 156.07803,3.7103437 158.0374,0.49038851 159.07253,0.093383369 C 160.15733,-0.32268751 160.58945,0.55956054 161.28726,4.6150731 C 161.64495,6.6938096 162.47416,9.0181138 163.30661,10.275395 C 166.28196,14.76915 169.6766,23.963119 170.34947,29.350168 C 170.56158,31.048448 171.09239,32.304761 171.93348,33.09925 C 173.26731,34.359168 174.16447,37.769092 174.8664,44.246714 C 175.32102,48.442199 175.4526,48.690133 178.11668,50.371026 L 180.18492,51.675976 L 179.89321,55.303463 C 179.61941,58.708303 179.70289,59.078872 181.25325,61.341065 C 183.20829,64.193738 183.30639,65.455783 181.80221,68.40272 C 180.44261,71.066398 180.92353,72.61609 183.52743,73.96193 C 186.55458,75.526528 191.68136,73.953162 196.88613,69.862236 C 198.01535,68.974674 198.51417,68.884664 199.73427,69.348307 C 202.52984,70.410636 203.70934,70.057868 205.46561,67.634164 C 207.63202,64.644475 211.79191,60.269908 212.46849,60.269908 C 215.41228,60.269908 215.95767,63.833895 213.48903,66.938704 C 210.10667,71.192661 209.33975,75.044942 210.54154,81.744199 C 210.86377,83.540431 210.70378,83.843073 208.15873,86.251487 C 205.95795,88.334122 205.43051,89.172844 205.43051,90.589831 C 205.46828,92.718038 207.72604,93.193174 208.54407,94.916659 C 208.81695,95.491593 207.60067,98.406018 206.2443,101.24959 C 203.06126,107.92274 202.81488,109.17763 204.29493,111.17848 C 206.75325,114.50183 206.76764,115.28752 204.45619,119.97909 C 199.6397,129.75513 199.53607,130.04012 199.53607,133.51009 C 199.53607,136.42216 199.72981,137.05805 201.07575,138.56365 C 203.1375,140.86997 205.50469,140.89925 208.40975,138.65439 C 209.57102,137.75704 211.51772,136.74356 212.73576,136.40225 C 214.36205,135.94653 215.61412,135.00514 217.44897,132.85852 C 219.41959,130.55307 220.61959,129.69047 223.12775,128.77642 C 224.87685,128.139 226.85521,127.04721 227.52411,126.35025 C 232.88129,120.76828 233.88343,120.32825 235.79651,122.71787 C 237.02952,124.3669 237.20709,126.37363 238.13203,128.18804 C 239.04762,129.68905 239.47855,129.93605 240.93511,129.79479 C 242.41824,129.65095 242.84269,129.24401 244.00408,126.85241 C 244.77647,125.26187 246.95329,122.50049 249.10287,120.38437 C 251.16593,118.35344 253.36463,116.0822 253.98888,115.33715 C 254.87747,114.27661 255.81944,113.90387 258.32771,113.62023 C 260.62879,113.36003 262.12357,112.82278 263.63267,111.71353 C 264.78831,110.86411 267.0276,109.67778 268.60887,109.07726 C 270.19014,108.47673 271.94069,107.49535 272.49897,106.8964 C 273.05726,106.29745 274.1705,105.67355 274.97283,105.50994 C 276.78542,105.14032 278.19153,104.36607 280.7279,102.34097 C 282.44312,100.9715 283.28452,100.72233 287.07053,100.46259 C 291.62076,100.15044 292.46386,99.771191 296.94679,96.020051 C 299.13884,94.185829 302.18774,94.019494 303.92721,95.63923 C 305.41725,97.026705 305.40698,97.176519 303.56953,100.85557 C 302.72778,102.54095 301.87,104.61939 301.66333,105.47431 C 301.45666,106.32924 300.57413,107.69485 299.70213,108.50901 C 298.83014,109.32315 297.77143,110.81516 297.34945,111.82456 C 296.92748,112.83398 295.70878,114.70477 294.64123,115.98189 C 293.24655,117.65035 292.79408,118.61621 293.0336,119.41344 C 294.01222,122.67059 293.96166,124.52486 292.84254,126.42066 C 291.47638,128.73494 291.43188,129.51896 292.56829,131.25246 C 293.47308,132.63264 293.41468,133.78992 292.22968,137.96103 C 291.27778,141.31166 291.78999,142.76242 294.6317,144.76442 C 297.64399,146.8866 297.73642,148.06749 295.18661,151.85456 C 294.15508,153.38661 293.27438,155.15817 293.22949,155.79135 C 293.1846,156.42454 293.09184,157.1836 293.02337,157.47818 C 292.95489,157.77274 292.87106,158.92294 292.83706,160.03416 C 292.78918,161.59914 292.2921,162.57846 290.63182,164.37871 C 288.22486,166.98861 288.17028,167.18999 289.28888,169.33367 C 290.80757,172.24416 290.44729,173.17561 285.65172,178.73682 C 283.56105,181.16129 283.02265,182.17835 283.02265,183.70324 C 283.02265,186.04288 284.31543,188.14035 286.19856,188.85594 C 288.63526,189.78191 292.74359,192.68658 293.06989,193.71413 C 293.55228,195.23322 292.87997,196.67202 291.18706,197.74355 C 289.29336,198.94217 284.15897,203.32832 281.24609,206.2358 C 280.05482,207.42486 277.39912,209.31339 275.34452,210.43254 C 270.48517,213.07947 270.20268,213.39108 269.68285,216.67774 C 269.16023,219.98218 268.07186,221.80092 265.89687,223.00445 C 262.88219,224.72434 260.55521,227.9158 256.80192,227.68991 C 251.83706,226.6619 252.03968,226.64571 250.41634,228.20017 C 248.06137,230.45521 246.95585,233.76371 247.75183,236.17429 C 248.68283,238.99382 250.54853,239.73174 257.6074,240.07237 C 261.90804,240.27991 265.04427,240.10767 268.85421,239.45473 C 271.68885,239.15253 274.76143,237.43506 275.87558,240.66957 C 276.36845,242.91242 277.64014,243.77783 280.66987,243.93217 C 283.65744,244.278 286.45854,245.0821 289.48878,245.12309 C 291.51267,245.09944 293.51819,245.40244 294.41885,245.86796 C 297.7468,246.83254 301.14988,246.6511 302.42068,250.03619 C 303.08041,251.79355 300.31743,253.48648 296.67058,254.16276 C 292.38955,254.95666 286.01902,257.96644 284.38629,259.96651 C 283.27293,261.33036 283.0417,261.37916 277.956,261.32352 C 275.05454,261.29178 271.69719,261.08684 270.49525,260.8681 C 267.55592,260.33316 266.28597,261.50532 266.67379,264.39525 C 267.0291,267.04294 266.19356,268.0778 263.39789,268.45259 C 260.6047,268.82703 259.98126,269.82646 260.44257,273.19028 C 260.75061,275.43652 260.64731,275.8217 259.5005,276.70293 C 258.78801,277.25041 257.31937,277.69801 256.21383,277.7046 C 251.81277,277.37462 248.84505,274.39105 244.56021,273.88878 C 244.11917,273.83708 242.97101,274.77021 241.99625,275.62562 C 239.8309,277.82314 238.65844,276.86361 236.12436,275.61604 C 230.19322,272.61929 224.56938,272.52007 223.84635,275.3994 C 223.50792,276.74708 223.84065,277.63443 225.89944,280.87474 L 227.05535,282.69402 L 225.82899,283.68655 C 224.1705,285.02883 223.67472,284.96702 218.75383,282.80454 C 216.40772,281.77355 214.24378,280.93 213.94509,280.93 C 213.64639,280.93 212.08263,280.30617 210.47008,279.54369 L 207.53814,278.15738 L 204.89781,279.10634 L 199.42318,278.65755 C 195.77759,276.85973 191.32009,276.40283 184.90028,277.16896 C 178.81996,277.89456 174.87257,277.14829 168.87323,274.13897 C 165.49436,272.44411 164.93441,272.32817 160.82615,272.47291 C 156.67731,272.61907 156.37074,272.71192 155.7497,274.01028 C 155.02077,275.53423 154.16019,283.47734 153.4333,295.39074 C 152.98575,302.72575 152.97137,307.23556 150.51525,317.88522 C 148.61389,326.12947 146.26136,335.53331 144.66752,338.56171 C 143.80163,340.20691 143.09316,341.89009 143.09316,342.30201 C 143.09316,343.1297 141.30651,343.24351 139.6101,342.5239 z\"/></svg>";
PShape leaf;
PShape meshShape;

void setup() {
  size(800, 800, FX2D);
  smooth();

  leaf = new PShapeSVG(parseXML(leafSVG));
  leaf = PGS_Transformation.scale(leaf, 1.4 * 1.6f);
  leaf = PGS_Transformation.rotateAroundCenter(leaf, -0.6f);
  leaf = PGS_Transformation.translateCentroidTo(leaf, width / 2f, height / 2f);
  leaf = PGS_Morphology.simplify(leaf, 1);

  meshShape = meshColorShape(leaf);
}

void draw() {
  background(0, 0, 40);
  shape(PGS_Morphology.fieldWarp(meshShape, 15, 0.5,frameCount/300f, false, 1337));
}

public void mouseMoved() {
  meshShape = meshColorShape(leaf);
}

PShape meshColorShape(PShape shape) {
  List<PVector> points = PGS_PointSet.poisson(0, 0, width, height, 45 + mouseX/33f, 10);
  IncrementalTin mesh = PGS_Triangulation.delaunayTriangulationMesh(shape, points, true, 1, true);
  
  String[] palette;
  if (mouseX > width/2) {
    meshShape = PGS_Meshing.urquhartFaces(mesh, true);
    palette = new String[] { "#563d7c", "#0096d8", "#f4e361", "#f24679" };
  } else {
    meshShape = PGS_Meshing.edgeCollapseQuadrangulation(mesh, true);
    palette = new String[] { "#07224f", "#ed361a", "#fc8405", "#f7c72a" };
  }

  PGS_Coloring.colorMesh(meshShape, ColoringAlgorithm.RLF, palette);
  PGS_Conversion.setAllStrokeColor(meshShape, 0, 1);
  PGS_Conversion.setAllStrokeToFillColor(meshShape);
  return meshShape;
}
