import processing.javafx.*;
import micycle.pgs.*;
import micycle.pgs.PGS_Coloring.ColoringAlgorithm;
import java.util.List;
import org.tinfour.common.IIncrementalTin;

String leafPolyline = "ai`tvCagmtmAojbVgpiw@oijDucv[b__AoflYp`xn@{xq|B|hbNsheT|pcTqyfP~es}@mxp\\|rcHm_tK|lbGchfp@sriAuys[u`oFutiYqdgKqypUe_rOieiPmn}Po{{IykyQywuCykrThd|@ulhYf_~GcchiAd~yg@q}yhA~hgq@iuovBf}oiBi~iSzniC{oqi@aosVs`cS{|nCklwUzhrC_qu|Bhn}{@ogshAl`ai@}yhfB`sza@}hsY_yRao}YalcKyezHygiPv_lCsv}WjjlOkyk\\vqgbAgocoA|fmc@ymcy@f~qVw~vs@~scFuatc@d_tGejtLjqgrBkusmA|dqaAslne@p{eJyrbQc_hAiazS}i|]yn{k@}saCe}fRjwbFm_~Qbclc@i~sg@zxmDcdaRozcFyg{Sgdpg@_qxt@`cbAa`sWdc{Kgq|TdrcS_viR~amXit_K||daCkz{\\hjmXoazI`tcH{fbOnfnPe|~u@jvdNihzLtewzAso~PxiccA_yhXddrWqkeOxc|Ou`tP~aeGc}lQg_q@capQ_lfRsisWukzsAiabg@||uMm}ei@xw_Me|eKn_jkDkpti@diauAqvxMvygX_h~a@~rkhAkczh@zk{dAwopSp_lo@{{hYtpdeBgounAxaccAoxrj@pudp@oogRj~pkC{xf[dkef@e~gSro~}@aben@`x|Li~iPfcIekxP_x}Mqrwe@unr|Aoa|jDk|hpBq`o~Dcrgq@gygbB{vmcAavggDstkn@adkaCyld\\_hfbBc~fNsbolAz|_G{}jQ`hga@kcuM??rfzg@wigGbvpKp`s{@x_xc@fprjBtpdeBp{lvFpw~cAbj_`Cjlz_BjrbiEtsbr@zj}}An{adAzdfqB|zbj@auoM~_tRkmuMbv`tBopbtBvw`j@mwad@fgleDoyytBhxh]wsh]tsnv@kcwwAj`dXmc~[nqygAtfIhnpUgtaF`orPwa}NdsuLuovX~mbScoux@fk~NityVphpz@mrnq@`ojx@ct~g@`d_|@m}sb@tfdWfxuKz{zDjdqSuuvInonTcrfAlw}_@d`mFdkdZvojQldnQ`wmTjswC|v_Y_daHto`_@}c`Vnzek@qeuj@ngu\\seld@hnfRsdcO`uqr@mpdXnvdgAkmfVbzk~@mw~gAz_mx@mbho@njrr@aezVnolR_sKxz|KbwpGtkePxgo|@fwnIzv}J|xoQpqcAtuzfAmw|Vt{gZoi{AjybUjikKbnvM|_eg@ndoGfh`GpsxMssH~cxU{}`Hn~miCgtusAf}pg@ye_SnxgxGq_}l@h}mr@zgyA~}jMl~xFzhkF~pdKuesAnxw[kzjIlzxJyh`uBhsfjAgvr`@tfi^gj~r@jng|@uk~~Avvwv@ezhNfxsPu`_IhlsTc}zLpabgA_tkrA|bmm@ogkxAnwsnAksrq@`mzY{|dUj|}Q}gcQvpqWwk~Kpji[mivExge]p}Rlje]xnoGtbyYz}aX~ekS~||Vhw_Kh~tx@bfnM~uao@}snArf`Sc{_Krglf@yqzc@l_bRqo{JrrgkAjkwCfmqrB`uqVtrlQ`hdG|w}h@tddq@vtaOldsJ~rcT`irDhgjl@luq@zk}pA_}}Jfen\\n_bBn|sv@r__Qdmxx@vahD`e}lB|myXlkbPbpfGtuw@bgzX{xtLj_vh@qixd@vkieAetjAdygVrw{Bfi|WhaiVrxnZvxwwC`xuaBn}b^nd}^~urUbyrn@xzuO``iTdazUld|Ptf|z@|voa@danVxhna@b`vRvlaSz~uu@dojf@t{o[tcsMjpdz@dkpSb_aVjgbQtrn@fqrN_yro@`inj@gtdGdhdT|oxAh|dUlruLbx_VxppWzwtVbr`oAvx{q@f_wLzfceApdqHrhnTvnwMtjjOb||m@vvq`@~jtNxqgy@lirJp}gUhziNnc`Nzqll@fxmQr}`vC|~xwAp_bqBlcuw@lcdNtzvS{psG||vSkpbUju~Ggd_fAz}kIk`{k@}wwFwgeq@vrkAohln@fg}Oy~leAj}wg@{xiPh_mAeo~lAgcgDizs}@zp|To__}@mqaDud|nBtgsF_vbtAh{\\o_gR`xiFet}hAv_jj@ycrLogcAyrjjAev}a@oqpeDg~l|@wvy`@gybFile]ka[szzS~kaEac~Jh{mQwpoF~wz]mhz@hoxa@hkjM`ob_Am~s@h`sT}p`LdzsMkipR|gfCosux@mjcFo}u]grqOayqu@e~mNg|bRrzdAej~c@vikSygdQrz}BsdmUoru@}mq_B_ss[y`zUadn@}lml@j}dGc`qhBaxdJk}aWdqjAykqTlaeHkhgQbykQavjDzxrQfagBlbrShhsKnz_WjxepDp{bgE~kre@rca]fpqGzt|MnmlF`||o@bs_Y~hsk@rrrp@dfii@rbkx@~hta@xb{x@xz}k@vrjShbbZwq{Lte_r@tvdIniwYdt~JdfcIljmn@rvcTlt~aAre}OxgzM`ynJzgd~Bz_toCtiim@ft_j@l|uUbvp[xrqKzxoXfet@x|{Saj`I`buMiwjTn|eCkks[knqDupsf@ypeScptp@_o}Jgnia@}wgB}lpVfzoCurzo@ngf[kiqXxagHg}hoB`nbXsbr|@v~xBy{fUlqlEytwLbl}P{nhCzu{Vl`uvAxypzEzhs\\fo{|@nkdGl|ja@sniAxoqRafqIhhdTy}~A~yq_@`hpRxjp~@fymNxdngAbzsPr{nb@tugFx~sn@qh{Eponq@hmkCthuVz}dtAtsrjEjrqHdmfdBpizJ|lbw@zroDrwlw@wqtDjvoa@c`eK~mvEgpkSyqpJc__t@exlq@wu_f@kzi[ekj~@as~]gnzkAglrn@icf_BgmejAukut@_y~s@umqPsihJyot_@u|mIomqPychIwszj@gb~d@uygrAwnvzA{ncRqecNyvhSoppDawyaAdhj@o_ux@mbk{AgvneAu`ib@qcsRkfmOkqpJqthUq{qLq|p{@urjL_ogQ_ysSertGco_\\f_Assqe@`~hNqsmc@loua@kxl_@nvur@spod@vdbwAsblo@lawQucoPlphL{ijH~{{R}jbLb|h|BwtbNrcay@aocUjzwFwueS_j{DuffO}hjP{beIc`sYswnBwbwzAy_lMck_s@crh\\mwjq@q_tz@y~neAgamAeqjT`rmLiuyjAsm\\_wxVkimG}jhQw`xUozxNoeiu@}naK{ppEc|_Jut}a@msplDkjtPurah@iafKe~oIoz_m@sd~Qq_dS{faLi~|Kw}dQ_gnGswdZ{_uYocfhDihdPy|pw@es}e@wuuw@i`fQcvqMwkfWy~fGqxyXks_BugdUbwsC{hsQdflK}|fN~ugS}jnTvprr@{mxf@nth_Ak}fHnbjZypjIh_vbAe}`Ptcbn@ch~]hzio@ovzQtd_k@i}uPzrxtAwkmH|{iUalqJjhlJs__NnlqAsypg@m`aM{giw@eugf@ejpZskrFcz~TtupJgasKdpcMal{Bp}pQftrC~|kq@{qjD~tfnAccuZ|wvxCehjKpghQwe`m@fvok@ibjOlzrUqhgc@|`alAi`mg@bundAc~oKvlpd@ke~c@pect@_qq[~qekAuzzIpmcPwhcvAv~zlA_wcWpwpf@_aja@|ko}AqkhSjmxTc|yWvxyKu~{Yfe~@kelUe_cEeoyImhlQe{mZiywpCrqwDg_vgA_jsCwx~g@f|_@_va|AcmsHajhR}hsj@up}k@m_mI_uzT}f_H_f}v@mtjJscyQ{hyd@oiuTqonKki_Nc}wf@eru~Aa{mPky}RcqhWulqIyc_~@geyDelsR}xtFklaK_v~MsedE}ajWsciBetz~@sq|F_lze@eyjh@qlglAil_Bap~Ul}yEs~weAeb~Eec_Sgp`i@mqvYslyK_vbLu|jFqjpPwprAmi|W|ocHqpmzB}c~AgpdY}hsHux~RobsP{xmRyifTabiLscaVigkEguqoBfbdFof~f@miyFakdMul|PwewBylzUt`o\\k|y_Bndob@wqqqCz{nS{euv@|mnh@om|yAbzxAaimSwobCgtoR";
PShape leaf;
PShape meshShape;

void setup() {
  size(1000, 1000, FX2D);
  smooth();

  leaf = PGS_Conversion.fromEncodedPolyline(leafPolyline);
  leaf = PGS_Transformation.translateCentroidTo(leaf, width / 2f, height / 2f);

  meshShape = meshColorShape(leaf);
}

void draw() {
  background(0, 0, 40);
  shape(PGS_Morphology.fieldWarp(meshShape, 15, 0.5, frameCount/300f, false, 1337));
}

public void mouseMoved() {
  meshShape = meshColorShape(leaf);
}

PShape meshColorShape(PShape shape) {
  List<PVector> points = PGS_PointSet.poisson(0, 0, width, height, 45 + mouseX/33f, 10);
  IIncrementalTin mesh = PGS_Triangulation.delaunayTriangulationMesh(shape, points, true, 1, true);

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
