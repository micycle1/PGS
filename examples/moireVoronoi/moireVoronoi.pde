import micycle.pgs.*;
import micycle.pgs.color.ColorUtils;
import processing.javafx.*;

String shape = "egj{iC}iwikDjnkl@gudHxnbi@b_iwCnof}BnlptHvxy|H`qiwQ|bwdAi`n_@pw{pDcqijDpxelEormsChc~eBim}qCnd~hBmnoF|s{iAyz{gC`otdEsfglCb}j`@`h~c@_~cE|ytyAnzfk@zmyXvvo_A{mtb@~wnfBsi`kBdnohCabxt@znfeBenxmBb{jxAsqmq@pl}b@dcaHtrw]fgspApxvaC_glZ`nqWzjwLnzjOr|sk@zt{XbdrGvsf_Fekl}BjyupJw`|n@rsdWxreU_q|Aza~^cnvaDtggpB}t`jAnsiuAgtsiBz`~|@_dpj@|hwzBgtzoGvn_~D_q`i@dti{An}cIlej~@`fnu@de_c@jb`}A|_lP|mwmBc}rfAnxvsAjihDbo{sChzgb@j~rm@xovv@ttjh@hwcR~gemD_cbGt`fjGn_v~@gmfLznajAiu~h@xmcmAl`kIpyfcAltioE~{nnCtz|w@rvb}Adq}zAfjgy@|fwkA|unrAzfbzC|s`oAywjs@v~s`Aog~Dzido@~doi@zr{r@~e}wAdxpw@za|XxgxdB|wpcA`tnv@nth]d}nxA|g~yJrc}|EchrHzu`V}pv_A|d_Qsd~iD`y`DevigBbxsw@ewinMljbXwaqfBb_iw@qy{zGqhnlB}ycfAmi{GgaeVf_wEq~sTdayu@rt`K|r|dCewec@|f~S{{{}Cggpj@_roqA|kn\\oyfcC}cv`@kud`E_{}C{|_fAp}g`@j`uJlsgeAxpxnFpnw~Fx_um@jtpwBbaylEvmrdDm~lNbqtw@x|hJ~mv\\nh~{C`ijaA`o~fFxxrkGi}hHzqpf@ovvVd|tCyphdBuwcc@k~`}AclkCgruqAvrsh@wk}xEtugf@mlaSd|nm@bojaCx`vwHiwqG|otvBjxh`Av|_wEaxn@pj|lBdzs}AfrfaFlssVle{rGilec@ufjD{ztoAclhhAs`~lEumzxBg~}hCwkcyB}qjjAgjrb@{apaDsw|_Dmiy_BudtCacb_AcqseBwhccB_iox@kz|Zcfa{Aijxd@m{|\\yxyt@l}uDeo}h@jwb[eop`BbzcpDyrfhApqrb@caoc@x}gvE}uhs@h_uKohj\\cvzo@}l{Rctv~Ca|yaBgoheChry@}swyCukttAco~]wiav@g`pgFgjdoBe{utAoscy@mgaxGc~o~@s|enAwesv@ouoKusil@ls}Qaz_tAjdivCuoh\\bwcdCgovz@ty`iB__hd@npolCsemr@n{{IqpisByfa~@{_ne@nrp[mym@l~cdCoqma@~}iiEqi`rApr|}Au|`wDhyhnJo}snB|goiBy}cq@rm~_Caxcq@nkne@mpxu@}prCeppi@i|luDpykAyxwdFgdez@wikfAcnfa@qrjjBku{v@edxg@kjr~@ouo_Cac_`Bc|eQ{qfb@uceXulgM_{v{Bqh`r@{eheBxcgCc|dgBgb|kAy~~tAhegFclzcDqh_M_ohr@iwxcA_loj@qmcxCpizCktb_@u`oXdx`cAiuzfH~vzcAy}{aDaaz]yiggDj}zdAk_bbDfxffBa`aaAfjlEq_puAun}Tgjiu@mhke@{lj^gq|k@_qbB{ssgDzkt{A_u}|Chl{gCieqtBg|iTa`jqHzshpCe`nxAezgE}p`Eej~m@~hmqBwdxxCdvx}@gv`rCn`ouEkrguCoplh@ytl}A|zlu@k}lvAiq~s@sj~qAfosO_cvr@rwsq@e{ya@bt|kD}xgl@tfmm@snk}Azu|mCihn_@|fzoA}xft@w`Ls}xv@ee{jBw`ex@rgjOo`|m@hxumHc_mlAfua[asze@`gwxD_`nqAv||}DoekjC`wi`GqohjA~t_kBuzmtAkx`h@{aboBitquFiecqLwzmyAa{txEap_eB}nftI~p~m@cuvc@";
PShape leaf;

void setup() {
  size(1000, 1000, FX2D);

  leaf = PGS_Conversion.fromEncodedPolyline(shape);
}

void draw() {
  background(0, 0, 40);

  shape(moireVoronoi(leaf));
}

PShape moireVoronoi(PShape shape) {
  var gridPoints = PGS_PointSet.squareGrid(0, 0, width, height, (10 + 181 / 10d));
  var grid2 = PGS_Conversion.fromPVector(PGS_PointSet.squareGrid(0, 0, width, height, (10 + 181 / 13d)));

  var c = PGS_ShapePredicates.centroid(shape);
  var gridRotate = PGS_Transformation.rotate(grid2, c, millis() / 2000f);
  gridPoints.addAll(PGS_Conversion.toPVector(gridRotate));

  var voronoi = PGS_Voronoi.innerVoronoi(gridPoints);
  var cells = PGS_Conversion.getChildren(voronoi);

  cells.forEach(cell -> {
    int id = Integer.parseInt(cell.getName());
    int col = ColorUtils.sinebow(id * 0.001 * 192);
    cell.setFill(col);
    cell.setStroke(col);
    cell.setStrokeWeight(0.5f);
  });

  return PGS_ShapeBoolean.intersectMesh(PGS_Conversion.flatten(cells), shape);
}
