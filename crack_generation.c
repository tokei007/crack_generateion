#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>

#define NUMBER_OF_CRACKFRONT_POINTS 10000
#define NUMBER_OF_CRACK_POINTS 100000
#define NUMBER_OF_PATCHS 1000000
#define NUMBER_OF_LAMINATION_POINTS 1000000
#define NUMBER_OF_TRIANGLE_VERTEXS 3
#define NUMBER_OF_OUTER_LAYERS 2
#define NUMBER_OF_INNER_LAYERS 2
#define DIMENSION 3
#define LAMINATIONSCALE 0.1
#define EPS 0.000001

#define CRACK_FRONT 2
#define CRACK_FACE 1
#define OUTER_CRACK 0

#define CRACK_SURFACE 1
#define OTHER 0

int nnodes;
int completed_nnodes;
int number_of_inner_nodes;
int number_of_innerest_nodes;
int npatch;
int minus_flag[NUMBER_OF_CRACKFRONT_POINTS];
int inout_flag[NUMBER_OF_CRACKFRONT_POINTS];
int crack_flag[NUMBER_OF_CRACK_POINTS];
int surface_flag[NUMBER_OF_CRACK_POINTS];
int new_nodes_flag[NUMBER_OF_LAMINATION_POINTS];
double nodes_coordinate[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double completed_nodes_coordinate[NUMBER_OF_CRACK_POINTS][DIMENSION];
double inner_nodes_coordinate[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double innerest_nodes_coordinate[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double patch_coord[NUMBER_OF_PATCHS][NUMBER_OF_TRIANGLE_VERTEXS][DIMENSION];
double patch_cg[NUMBER_OF_PATCHS][DIMENSION];
double init_univec_normal[DIMENSION];
double univec_normal[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double univec_propa[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double univec_tangent[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double innerest_univec_normal[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double innerest_univec_propa[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double innerest_univec_tangent[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double p_to_p_vector[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double innerest_nodes_p_to_p_vector[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
double layer_size;
double scale_factor;
double lamination_weight[NUMBER_OF_CRACKFRONT_POINTS];
double layer_length[NUMBER_OF_CRACKFRONT_POINTS];
double layer_length_sum;
char patch_file_name[256];
double temp_crackface_new_nodes_coordinate[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
int temp_count;


double Distance(double a[], double b[])
{
  double dist;
  //printf("##a[0,1,2] = %lf %lf %lf, b[0,1,2] = %lf %lf %lf\n", a[0], a[1], a[2], b[0], b[1], b[2]);
  dist = sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]));
  //printf("##dist = %lf\n", dist);
  return dist;
}

double InnerProduct(double a[], double b[])
{
  double product;
  product = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  return product;
}

/*
 * ベクトルa,bの外積ベクトルcを算出
 * (a×b=c)
 */
void CrossProduct(double a[], double b[], double c[])
{
  //printf("##a[0,1,2] = %lf %lf %lf, b[0,1,2] = %lf %lf %lf\n", a[0], a[1], a[2], b[0], b[1], b[2]);
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

double GetVectorLength(double a[])
{
  double length;
  length = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  return length;
}

void GetMaxandMin(int start, int total, double value[], int max, int min)
{
  int i;
  double max_value;
  double min_value;
  max_value = value[start];
  max = start;
  min_value = value[start];
  min = start;
  for(i = 0; i < total; i++){
    if(max_value < value[i]){
      max_value = value[i];
      max = i;
      if(min_value > value[i]){
        min_value = value[i];
        min = i;
      }
    }
  }
}

int InverseMatrix_3D( double M[3][3], double zero ){
  int i, j;
  double a[3][3];
  double det = M[0][0]*M[1][1]*M[2][2] +M[0][1]*M[1][2]*M[2][0] +M[0][2]*M[1][0]*M[2][1]
    -M[0][0]*M[1][2]*M[2][1] -M[0][2]*M[1][1]*M[2][0] -M[0][1]*M[1][0]*M[2][2];

  //	printf("det = %20.15e\n",det);
  if(fabs(det) < zero) return(1); // matrix is singular 

  for( i = 0; i < 3; i++ ){
    for( j = 0; j < 3; j++ )	a[i][j] = M[i][j];
  }
  M[0][0] = (a[1][1]*a[2][2]-a[1][2]*a[2][1])/det; M[0][1] = (a[0][2]*a[2][1]-a[0][1]*a[2][2])/det; M[0][2] = (a[0][1]*a[1][2]-a[0][2]*a[1][1])/det;
  M[1][0] = (a[1][2]*a[2][0]-a[1][0]*a[2][2])/det; M[1][1] = (a[0][0]*a[2][2]-a[0][2]*a[2][0])/det; M[1][2] = (a[0][2]*a[1][0]-a[0][0]*a[1][2])/det;
  M[2][0] = (a[1][0]*a[2][1]-a[1][1]*a[2][0])/det; M[2][1] = (a[0][1]*a[2][0]-a[0][0]*a[2][1])/det; M[2][2] = (a[0][0]*a[1][1]-a[0][1]*a[1][0])/det;
  return (0);
}
void GetAverageVector(double a[], double b[], double c[])
{
  int i;
  for(i = 0; i < DIMENSION; i++){
    c[i] = (a[i] + b[i])/2;
  }
}

void GetUnitVector(double a[])
{
  double temp_amount;
  temp_amount = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

  a[0] = a[0]/temp_amount;
  a[1] = a[1]/temp_amount;
  a[2] = a[2]/temp_amount;
}

void ReadNodes(const char *filename)
{
  int i,nodeid[NUMBER_OF_CRACKFRONT_POINTS];
  FILE *fp;

  fp = fopen (filename,"r");
  assert(fp!=NULL);

  fscanf(fp,"%d \n",&nnodes);

  for(i=0;i<nnodes;i++)
  {
    fscanf(fp,"%d %le %le %le\n",&nodeid[i],&nodes_coordinate[i][0],&nodes_coordinate[i][1],&nodes_coordinate[i][2]);
    printf("x=%le y=%le z=%le\n",nodes_coordinate[i][0],nodes_coordinate[i][1],nodes_coordinate[i][2]);
  }
}


void ReadUnivector(const char *filename)
{
  FILE *fp;

  fp = fopen (filename,"r");
  assert(fp!=NULL);

  fscanf(fp,"%le %le %le\n",
      &init_univec_normal[0],
      &init_univec_normal[1],
      &init_univec_normal[2]);
  fclose(fp);
}

void ReadPatch(const char *fileName)
{
  //static int npatch;
  //static double patch_coord[MAX_PATCH][3][3];

  FILE *fp;
  int ipatch, ip, idir;
  int idummy;

  fp = fopen (fileName, "r");

  fscanf(fp,"%d", &npatch);

  for(ipatch=0; ipatch < npatch; ipatch++)
  {
    fscanf(fp,"%d",&idummy);
    for(ip=0; ip < 3; ip++)
      for(idir=0; idir < 3; idir++)
        fscanf(fp,"%lf",&patch_coord[ipatch][ip][idir]);
  }
  fclose(fp);
}

#define INNEREST "___temp_innerest.node"
#define INNERESTINOUT "___temp_innerest.node_inout"

void WriteNodes(int total, double nodes_coordinate[][DIMENSION])
{
  const char filename[200] = "___temp_innerest.node";
  FILE *fp;
  int inode,idir;

  fp = fopen(filename, "w");

  fprintf(fp, "%d\n", total);

  for(inode = 0; inode < total; inode++){
    fprintf(fp, "%d %lf %lf %lf\n", inode, nodes_coordinate[inode][0],nodes_coordinate[inode][1],nodes_coordinate[inode][2]);
  }
  fclose(fp);
}

void ReadNodesInOut(int total)
{
  const char *filename = "INNERESTINOUT";
  FILE *fp;
  int temp_total;
  int i;
  int dummy;
  fp = fopen(filename, "r");

  fscanf(fp, "%d", &temp_total);
  assert(total==temp_total);

  for(i = 0; i < temp_total; i++){
    fscanf(fp, "%d %d", &dummy, &inout_flag[i]);
  }
}

void WriteCompleteNodes(const char *filename)
{
  FILE *fp;
  int inode,idir;

  fp = fopen(filename, "w");

  fprintf(fp, "%d\n", completed_nnodes);

  for(inode = 0; inode < completed_nnodes; inode++){
    fprintf(fp, "%d %lf %lf %lf\n", inode, completed_nodes_coordinate[inode][0],completed_nodes_coordinate[inode][1],completed_nodes_coordinate[inode][2]);
  }
  fclose(fp);
}

void WriteCrackFlag(const char *filename)
{
  FILE *fp;
  int inode;

  fp = fopen(filename, "w");

  fprintf(fp, "%d\n", completed_nnodes);

  for(inode = 0; inode < completed_nnodes; inode++){
    fprintf(fp, "%d %d\n", inode, crack_flag[inode]);
  }
  fclose(fp);
}

void WriteSurfaceFlag(const char *filename)
{
  FILE *fp;
  int inode;

  fp = fopen(filename, "w");

  fprintf(fp, "%d\n", completed_nnodes);

  for(inode = 0; inode < completed_nnodes; inode++){
    fprintf(fp, "%d %d\n", inode, surface_flag[inode]);
  }
  fclose(fp);
}

void AddNewNode(double new_nodes_coordinate[], int temp_crack_flag, int temp_surface_flag)
{
  int i;
  for(i = 0; i < DIMENSION; i++){
    completed_nodes_coordinate[nnodes][i] = new_nodes_coordinate[i];
  }
  crack_flag[nnodes] = temp_crack_flag;
  surface_flag[nnodes] = temp_surface_flag;
  completed_nnodes++;
}

void TempAddNewNode(double new_nodes_coordinate[])
{
  int i;
  for(i = 0; i < DIMENSION; i++){
  temp_crackface_new_nodes_coordinate[temp_count][i] = new_nodes_coordinate[i];
  }
  new_nodes_flag[temp_count] = 1;
  temp_count++;
}

void ClearNumberOfInnerestNodes()
{
  number_of_innerest_nodes = 0;
}

void ClearInOutFlag()
{
  int i;
  for(i = 0; i < NUMBER_OF_CRACKFRONT_POINTS; i++){
    inout_flag[i] = 1;
  }
}

void ClearNewNodeFlag()
{
  int i;
  for (i = 0; i <NUMBER_OF_LAMINATION_POINTS; i++)
    new_nodes_flag[i] = 0;
}

void AddInnerNode(double new_nodes_coordinate[])
{
  int i;
  for(i = 0; i < DIMENSION; i++){
    inner_nodes_coordinate[number_of_inner_nodes][i] = new_nodes_coordinate[i];
  }
  number_of_inner_nodes++;
}

void ResisterCrackFrontNodes(int total, double temp_nodes_coordinate[][DIMENSION])
{
  int i;
  for(i = 0; i < total; i++){
    if(i == 0 || i == total-1){
      AddNewNode(temp_nodes_coordinate[i], CRACK_FRONT, CRACK_SURFACE);
    } else {
      AddNewNode(temp_nodes_coordinate[i], CRACK_FRONT, OTHER);
    }
  }
}

//patch_cgとして表面パッチの重心を登録
void CompPatchCG(void)
{
  int iPatch, inode, idir;

  for(iPatch=0; iPatch < npatch; iPatch++)
    for(idir=0; idir<3; idir++)
      patch_cg[iPatch][idir] 
        = (patch_coord[iPatch][0][idir]
            + patch_coord[iPatch][1][idir] + patch_coord[iPatch][2][idir])/3.0;

}

//ThreshValとしてすべてのパッチで最も長い辺の＊５をreturn
double CompThreshDist(void)
{
  // find the smappset patch;

  double min = 1000000000000.0;
  double max = 0.0;
  int iPatch;
  int idir,i0, i1;
  double TempVec[3], Abs_value;
  double ThreshVal;

  for(iPatch=0; iPatch < npatch; iPatch++)
  {
    for(i0=0; i0<3; i0++)
    {
      if(i0<2) i1=i0+1;
      if(i0==2) i1 = 0;
      for(idir=0; idir<3; idir++)
        TempVec[idir] 
          = patch_coord[iPatch][i1][idir] - patch_coord[iPatch][i0][idir];

      Abs_value = sqrt(TempVec[0]*TempVec[0] + TempVec[1]*TempVec[1]
          + TempVec[2]*TempVec[2]);
      if(Abs_value < min) min = Abs_value;
      if(Abs_value > max) max = Abs_value;	  
    }
  }

  //  printf("Patch size: min and max: %15.5e and %15.5e\n",min,max);
  ThreshVal = max*5;
  return(ThreshVal);
}

//ある点から出ているベクトルと三角形パッチが交わるときの
//vecA,vecB,-vecTの係数を引数として出力
int CompCrossingPt_TriangleAndLine(double vecA[3], double vecB[3], 
    double vecT[3], double xxT[3], double xx0[3],
    double *pp, double *qq, double *mm, double zero)
{
  double matA[3][3], rhsv[3], solvec[3];
  int ii, jj;

  for(ii=0; ii<3; ii++){
    matA[ii][0] = vecA[ii];
    matA[ii][1] = vecB[ii];
    matA[ii][2] = -vecT[ii];

    rhsv[ii] = xx0[ii] - xxT[ii];
  }

  if(InverseMatrix_3D(matA, zero) == 1) return(1); //matrix is singular

  for(ii=0; ii<3; ii++) solvec[ii]=0.0;
  for(ii=0; ii<3; ii++)
    for(jj=0; jj<3; jj++)
      solvec[ii] += matA[ii][jj] * rhsv[jj];

  *pp = solvec[0];
  *qq = solvec[1];
  *mm = solvec[2];

  return(0);
}

void GetLayerLength(int total, double nodes_coordinate[][DIMENSION])
{
  int i, j;
  double value;
  double sum;
  double temp_vector[DIMENSION];
  for(i = 0; i < total - 1; i++){
    for(j = 0; j < DIMENSION; j++){
      temp_vector[j] = nodes_coordinate[i+1][j] - nodes_coordinate[i][j];
    }
    layer_length[i] = sqrt(temp_vector[0]*temp_vector[0] + temp_vector[1] * temp_vector[1] + temp_vector[2] * temp_vector[2]);
    layer_length_sum += layer_length[i];
  }
}

void RearrangeLayer(int total, double nodes_coordinate[][DIMENSION])
{
  int i, j;
  int number_of_division;
  int count = 1;
  int nodes_count = 0;
  double temp_nodal_distance;
  double vector_coefficient = 0.0;
  double temp_nodes_coordinate[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];
  number_of_division = layer_length_sum/layer_size + 0.5;
  temp_nodal_distance = layer_length_sum/number_of_division;
  for(j = 0; j < DIMENSION; j++){
    temp_nodes_coordinate[0][j] = nodes_coordinate[0][j];
  }

  while(1){
    vector_coefficient = vector_coefficient - temp_nodal_distance;
    if(vector_coefficient >= 0){
      for(j = 0; j < DIMENSION; j++){
        temp_nodes_coordinate[count][j] = nodes_coordinate[nodes_count][j] + vector_coefficient * (nodes_coordinate[nodes_count+1][j] - nodes_coordinate[nodes_count][j]);
        count++;
      }
    } else {
      vector_coefficient = vector_coefficient + layer_length[nodes_count];
      nodes_count++;
      if(nodes_count == total) break;
    }
  }
  for(j = 0; j < DIMENSION; j++){
    temp_nodes_coordinate[count][j] = nodes_coordinate[total][j];
  }
  for(i = 0; i < count; i++){
    for(j = 0; j < DIMENSION; j++){
      nodes_coordinate[i][j] = temp_nodes_coordinate[i][j];
    }
  }
  total = count;
}

//き裂前縁メッシュサイズ及び層間距離読み込み
void ReadCrackParam(const char *fileName)
{
  FILE *fp;
  double d1,d2,d3;

  fp = fopen (fileName, "r");

  if (fp == NULL) {
    fprintf (stderr, " file %s not found \n", fileName);
    exit(1);
  }

  fscanf(fp,"%lf %lf",&d1,&d2);
  fscanf(fp,"%lf",&d3);

  layer_size = d3;
}

void ReadVector(const char *filename)
{
}

void GetPointtoPointVector(int temp_nnodes, double temp_nodes_coordinate[][DIMENSION], double temp_p_to_p_vector[][DIMENSION])
{
  int i, j;
  for(i = 0; i < temp_nnodes-1; i++){
    for(j = 0; j < DIMENSION; j++){
      temp_p_to_p_vector[i][j] = temp_nodes_coordinate[i+1][j] - temp_nodes_coordinate[i][j];
      GetUnitVector(temp_p_to_p_vector[i]);
    }
  }
}

void SetScaleFactor()
{
  scale_factor = layer_size * 2;
}

void ClearMinusFlag()
{ int i;
  for(i = 0; i < NUMBER_OF_CRACKFRONT_POINTS; i++){
    minus_flag[i] = 0;
  }
}

void GetUnivectoratAllPoint(int temp_nnodes, double temp_p_to_p_vector[][DIMENSION], double temp_univec_normal[][DIMENSION], double temp_univec_propa[][DIMENSION], double temp_univec_tangent[][DIMENSION])
{
  int i, j;
  double product;
  for(j = 0; j < DIMENSION; j++){
    temp_univec_tangent[0][j] = temp_p_to_p_vector[0][j];
    temp_univec_normal[0][j] = init_univec_normal[j];
  }
  GetUnitVector(temp_univec_normal[0]);
  CrossProduct(temp_univec_normal[0], temp_univec_tangent[0], temp_univec_propa[0]);
  GetUnitVector(temp_univec_propa[0]);


  for(i = 1; i < temp_nnodes - 1; i++){
    CrossProduct(temp_p_to_p_vector[i], temp_p_to_p_vector[i-1], temp_univec_normal[i]);
    product = InnerProduct(temp_univec_normal[i], temp_univec_normal[i-1]);
    if(product < 0){
      for(j = 0; j < DIMENSION; j++){
        temp_univec_normal[i][j] = -temp_univec_normal[i][j];
      }
      minus_flag[i] = 1;
    }
    GetAverageVector(temp_p_to_p_vector[i-1], temp_p_to_p_vector[i], temp_univec_tangent[i]);
    GetUnitVector(temp_univec_normal[i]);
    CrossProduct(temp_univec_normal[i], temp_univec_tangent[i], temp_univec_propa[i]);
    GetUnitVector(temp_univec_propa[i]);
  }

  for(j = 0; j < DIMENSION; j++){
    temp_univec_tangent[temp_nnodes-1][j] = temp_p_to_p_vector[temp_nnodes-2][j];
    temp_univec_normal[temp_nnodes-1][j] = temp_univec_normal[temp_nnodes-2][j];
  }
  CrossProduct(temp_univec_normal[temp_nnodes-1], temp_univec_tangent[temp_nnodes-1], temp_univec_propa[temp_nnodes-1]);
  GetUnitVector(temp_univec_propa[temp_nnodes-1]);
}

void SurfaceCorrectionAdvVec(double AdvVec[3], double Normal[3])
{
  double InnerProd;

  InnerProd = AdvVec[0]*Normal[0]+AdvVec[1]*Normal[1]+AdvVec[2]*Normal[2];
  AdvVec[0] = AdvVec[0] - InnerProd * Normal[0];
  AdvVec[1] = AdvVec[1] - InnerProd * Normal[1];
  AdvVec[2] = AdvVec[2] - InnerProd * Normal[2];

  //  printf("Check corrected vec inner product=%10.3e\n",AdvVec[0]*Normal[0]+
  //	 + AdvVec[1]*Normal[1]+AdvVec[2]*Normal[2]);
}

void ComputeNormalNormalAndMovingVec(double CoordVec[3], double TanVec[3],
    double Normal[3], double MovVec[3], 
    double DistVec)
{

  int idir, iPatch;
  double pt_coord[3], pt_CG[3], dist1;
  double vecA[3], vecB[3], vecT[3], xx0[3], xxT[3];
  double pp,qq,mm;
  double zero = 0.000000000000001;
  double epsilon = 0.0000001;
  int if_singular;
  double Min_mm = 100000000000.0;
  double pp_keep, qq_keep, mm_keep;
  int iPatch_keep, iCount = 0;
  double ThreshDist;

  //patch_cgとして表面パッチの重心を登録
  CompPatchCG();
  //すべてのパッチで最も長い辺の＊５をreturnする関数
  ThreshDist = CompThreshDist();

  for(idir=0; idir<3; idir++){
    pt_coord[idir] = CoordVec[idir];
    xx0[idir] = CoordVec[idir]; 
    vecT[idir] = TanVec[idir];
  }

  for(iPatch=0; iPatch < npatch; iPatch++){
    pt_CG[0] = patch_cg[iPatch][0];
    pt_CG[1] = patch_cg[iPatch][1];
    pt_CG[2] = patch_cg[iPatch][2];

    //パッチの重心から判定したい点までの距離を算出
    dist1 = sqrt((pt_CG[0]-pt_coord[0])*(pt_CG[0]-pt_coord[0])
        +(pt_CG[1]-pt_coord[1])*(pt_CG[1]-pt_coord[1])
        +(pt_CG[2]-pt_coord[2])*(pt_CG[2]-pt_coord[2]));

    //もしThreshDistより小さければ
    if(dist1 < ThreshDist){
      for(idir=0; idir<3; idir++){
        vecA[idir] = patch_coord[iPatch][1][idir] 
          - patch_coord[iPatch][0][idir];
        vecB[idir] = patch_coord[iPatch][2][idir] 
          - patch_coord[iPatch][0][idir];
        xxT[idir] = patch_coord[iPatch][0][idir];
      }

      //き裂前縁の端の点から出ているベクトルと三角形パッチが交わるときの
      //vecA,vecB,-vecTの係数を引数として出力
      if_singular 
        = CompCrossingPt_TriangleAndLine(vecA, vecB, vecT, xxT, xx0,
            &pp, &qq, &mm, zero);
      //fprintf(fp,"Patch# = %d singul.?=%d  pp = %10.3e qq = %10.3e mm = %10.3e\n",	 
      //	  iPatch, if_singular, pp, qq, mm);

      if(if_singular == 0){
        ////き裂前縁の端の点から出ているベクトルと三角形パッチが交わるときの
        //vecA,vecB,-vecTの係数の係数をみて、三角形が成す平面と
        //ベクトルの交点がパッチの面内にあり、点からパッチまでの距離が最も近いところをiPatch_keepとする
        if(Min_mm > fabs(mm) && (pp+qq) < (1.0+epsilon) && 
            pp > (-epsilon) && qq > (-epsilon)){
          Min_mm = fabs(mm);
          pp_keep = pp;
          qq_keep = qq;
          mm_keep = mm;
          iPatch_keep = iPatch;
          iCount++;
        }
      }
    }
  }

  //  printf("Patch# = %d pp = %10.3e qq = %10.3e mm = %10.3e Count=%d\n",
  //	 iPatch_keep, pp_keep, qq_keep, mm_keep, iCount);

  if(iCount == 0){
    printf("Error counter is ZERO: Serching in MovingCoordinate\n");
    exit(1);
  }   
  for(idir=0; idir<3; idir++){
    vecA[idir] = patch_coord[iPatch_keep][1][idir] 
      - patch_coord[iPatch_keep][0][idir];
    vecB[idir] = patch_coord[iPatch_keep][2][idir] 
      - patch_coord[iPatch_keep][0][idir];
    MovVec[idir] = mm_keep * vecT[idir];
  }
  //パッチが成す面の法線ベクトル（平面の方程式の導出に使用）を算出
  CrossProduct(vecA, vecB, Normal);
  //法線ベクトルを単位ベクトルに
  GetUnitVector(Normal);
  //DistVecに点からパッチまでの距離を登録
  DistVec = sqrt(MovVec[0]*MovVec[0]+MovVec[1]*MovVec[1]+MovVec[2]*MovVec[2]);

  //  printf("Normal Vec  %10.3e %10.3e   %10.3e\n",Normal[0],Normal[1],Normal[2]);
  //  printf("Moving Distance =  %10.3e\n",DistVec);
}

void VectorCorrectiontoSurface()
{

  int i,j;
  double CoordVecStart[3], CoordVecEnd[3]; 
  double AdvCoordVecStart1[3], AdvCoordVecEnd1[3];
  double AdvCoordVecStart2[3], AdvCoordVecEnd2[3];
  double TanVecStart[3], TanVecEnd[3];
  double AdvVecStart[3], AdvVecEnd[3];
  double NormalStart[3], NormalEnd[3];
  double MovVecStart[3], MovVecEnd[3];
  double AdvVecStart1[3], AdvVecEnd1[3];
  double AdvVecStart2[3], AdvVecEnd2[3];
  double weight1=0.5, weight2=0.5;
  double DistVecStart, DistVecEnd;

  //き裂前縁端の点２点をStartとEndとして登録
  //それぞれ、隣の点から端の点へと向かうベクトルをTanVecとして登録
  //進展ベクトルにscale_factorを掛けたものをAdvVecとして登録
  for(i=0; i<3; i++){
    CoordVecStart[i] = nodes_coordinate[0][i];
    CoordVecEnd[i] =  nodes_coordinate[nnodes-1][i];
    TanVecStart[i] = -univec_tangent[0][i];
    TanVecEnd[i] = univec_tangent[nnodes-1][i];
    AdvVecStart[i] = univec_propa[0][i]*scale_factor;
    AdvVecEnd[i] = univec_propa[nnodes-1][i]*scale_factor;
  }  

  //TanVecを単位ベクトルに変換
  GetUnitVector(TanVecStart);
  GetUnitVector(TanVecEnd);

  //き裂の端の点について、Tanvecの方向に存在する最も近いパッチが成す
  //法線ベクトルと、端の点からパッチを結ぶベクトルと、その長さを導出
  //(ただし前ステップのき裂前縁点はそもそも表面に存在するので
  //MovVecとDistVecはほぼ無意味)
  ComputeNormalNormalAndMovingVec(CoordVecStart,TanVecStart,NormalStart,
      MovVecStart, DistVecStart);
  ComputeNormalNormalAndMovingVec(CoordVecEnd,TanVecEnd,NormalEnd,
      MovVecEnd, DistVecEnd);

  //進展ベクトルを法線ベクトルに対して垂直になるように補正
  //(進展ベクトルの法線に垂直な方向の成分を取り出す感じ)
  SurfaceCorrectionAdvVec(AdvVecStart, NormalStart);
  SurfaceCorrectionAdvVec(AdvVecEnd, NormalEnd);

  for(i=0;i<3;i++){
    //き裂前縁の端の点をAdvVecだけ動かした点
    AdvCoordVecStart1[i] = CoordVecStart[i] + AdvVecStart[i];
    AdvCoordVecEnd1[i] = CoordVecEnd[i] + AdvVecEnd[i];
    TanVecStart[i] = NormalStart[i];
    TanVecEnd[i] = NormalEnd[i];
  }

  //き裂の端の点からAdvVecだけ進んだ点について、Tanvecの方向に存在する最も近いパッチが成す
  //法線ベクトルと、端の点からパッチを結ぶベクトルと、その長さを導出
  ComputeNormalNormalAndMovingVec(AdvCoordVecStart1,TanVecStart,NormalStart,
      MovVecStart, DistVecStart);
  ComputeNormalNormalAndMovingVec(AdvCoordVecEnd1,TanVecEnd,NormalEnd,
      MovVecEnd, DistVecEnd);

  //き裂の端の点からAdvVecだけ進んだ点を最も近いパッチ面上へ移動
  for(i=0; i<3; i++){
    AdvCoordVecStart2[i] = AdvCoordVecStart1[i] + MovVecStart[i];
    AdvCoordVecEnd2[i] =  AdvCoordVecEnd1[i] + MovVecEnd[i];
  } 

  //き裂端の点から垂直補正した後のAdvVecだけ進んだ点までのベクトルと
  //き裂端の点から垂直補正した後のAdvVecだけ進んだ点からパッチ上へ移動させた点へのベクトルを登録
  for(i=0; i<3; i++){
    AdvVecStart1[i] = AdvCoordVecStart1[i] - CoordVecStart[i];
    AdvVecEnd1[i] = AdvCoordVecEnd1[i] - CoordVecEnd[i];
    AdvVecStart2[i] = AdvCoordVecStart2[i] - CoordVecStart[i];
    AdvVecEnd2[i] = AdvCoordVecEnd2[i] - CoordVecEnd[i];
  }

  //重みを半々として端の点における進展ベクトルを登録
  for(i=0; i<3; i++){
    univec_propa[0][i]= (AdvVecStart1[i]*weight1+AdvVecStart2[i]*weight2)
      /scale_factor;
    univec_propa[nnodes-1][i] =(AdvVecEnd1[i]*weight1+AdvVecEnd2[i]*weight2)
      /scale_factor;
  }
}

void MakeInOutLayer()
{
  int i,j,k;
  double coord_of_outer_layer[DIMENSION];
  double coord_of_inner_layer[DIMENSION];
  for(i = 0; i < NUMBER_OF_OUTER_LAYERS; i++){
    for(j = 0; j < nnodes; j++){
      for(k = 0; k < DIMENSION; k++){
        coord_of_outer_layer[k] = nodes_coordinate[j][k] + (i+1) * univec_propa[j][k];
      }
      if(j == 0 || j == nnodes-1){
        AddNewNode(coord_of_outer_layer, OUTER_CRACK, CRACK_SURFACE);
      } else {
        AddNewNode(coord_of_outer_layer, OUTER_CRACK, OTHER);
      }
    }
  }

  ClearNumberOfInnerestNodes();

  for(i = 0; i < NUMBER_OF_INNER_LAYERS; i++){
    for(j = 0; j < nnodes - 1; j++){
      for(k = 0; k < DIMENSION; k++){
        coord_of_inner_layer[k] = nodes_coordinate[j][k] - (i+1) * univec_propa[j][k];
      }
      if(j == 0 || j == nnodes-1){
        AddNewNode(coord_of_inner_layer, CRACK_FACE, CRACK_SURFACE);
      } else {
        AddNewNode(coord_of_inner_layer, CRACK_FACE, OTHER);
      }
      if(i == NUMBER_OF_INNER_LAYERS-1){
        AddInnerNode(coord_of_inner_layer);
      }
    }
  }
}

void DecisionOfInOut(double temp_nnodes, double temp_nodes_coordinate[][DIMENSION])
{
  char buf[256];
  WriteNodes(temp_nnodes, temp_nodes_coordinate);
  sprintf(buf, "~/ADVENTURE/advauto_etc_new/autoage/hetare/gm3d/advautoage_h_gm3d_tri_inout2 %s %s %s", patch_file_name, INNEREST, INNERESTINOUT);
  ReadNodesInOut(temp_nnodes);
}

void CalculateLaminationWeight(int temp_nnodes, double temp_p_to_p_vector[][DIMENSION], double temp_lamination_weight[])
{
  int i;
  int max, min;
  double temp_inner_product;
  double temp_cross_product[DIMENSION];
  double temp_cross_product_length;
  double temp_value[NUMBER_OF_CRACKFRONT_POINTS];
  double temp_temp_value[NUMBER_OF_CRACKFRONT_POINTS];
  for(i = 1; i < temp_nnodes - 1; i++){
    CrossProduct(temp_p_to_p_vector[i], temp_p_to_p_vector[i-1], temp_cross_product);
    temp_inner_product = InnerProduct(temp_p_to_p_vector[i], temp_p_to_p_vector[i-1]);
    temp_cross_product_length = GetVectorLength(temp_cross_product);
    if(temp_inner_product < 0){
      temp_value[i] = temp_cross_product_length;
    } else {
      temp_value[i] = 2 - temp_cross_product_length;
    }
    if(minus_flag[i] == 1){
      temp_value[i] = -temp_value[i];
    }
  }
  GetMaxandMin(1, temp_nnodes-1, temp_value, max, min);
  for(i = 1; i < temp_nnodes - 1; i++){
    temp_temp_value[i] = temp_value[i] - temp_value[min] + (temp_value[max] - temp_value[min])/10;
    temp_lamination_weight[i] = temp_temp_value[i]/(temp_value[max] - temp_value[min]);
  }
}

void LaminationLayer()
{
  int i,j;
  int inside_count;
  double lamination_scale_factor = layer_size * LAMINATIONSCALE;
  double temp_innerest_nodes_coordinate[NUMBER_OF_CRACKFRONT_POINTS][DIMENSION];

  number_of_innerest_nodes = number_of_inner_nodes;
  for(i = 0; i < number_of_inner_nodes; i++){
  innerest_nodes_coordinate[i][j] = inner_nodes_coordinate[i][j];
  }
  temp_count = 0;
  ClearNewNodeFlag();

  while(1){
    ClearMinusFlag();
    GetPointtoPointVector(number_of_innerest_nodes, innerest_nodes_coordinate, innerest_nodes_p_to_p_vector);
    GetUnivectoratAllPoint(number_of_innerest_nodes, innerest_nodes_p_to_p_vector, innerest_univec_normal, innerest_univec_propa, innerest_univec_tangent);
    CalculateLaminationWeight(number_of_innerest_nodes, innerest_nodes_p_to_p_vector, lamination_weight);
    for(i = 1; i < number_of_innerest_nodes - 1; i++){
      if(inout_flag[i] == 1){
        for(j = 0; j < DIMENSION; j++){
          innerest_nodes_coordinate[i][j] = innerest_nodes_coordinate[i][j] - lamination_scale_factor * lamination_weight[i] * innerest_univec_propa[i][j];
        }
        inside_count++;
      }
    }
    GetLayerLength(number_of_innerest_nodes, innerest_nodes_coordinate);
    for(i = 0; i < number_of_innerest_nodes - 1; i++){
      if(layer_size * 0.5 > layer_length[i] || layer_size * 1.1 <layer_length[i]){
        RearrangeLayer(number_of_innerest_nodes, innerest_nodes_coordinate);
        break;
      }
    }
    for(i = 1; i < number_of_innerest_nodes - 1; i++){
    TempAddNewNode(innerest_nodes_coordinate[i]);
    }
    if(inside_count == 0) break;
    ClearInOutFlag();
    DecisionOfInOut(number_of_innerest_nodes, innerest_nodes_coordinate);
  }
  for(i = 1; i < number_of_innerest_nodes-1; i++){
    AddNewNode(innerest_nodes_coordinate[i], CRACK_FACE, CRACK_SURFACE);
  }
}

void LaminationPointsToNodes()
{
  int i, j;
  double dist;
  for(i = 0; i < number_of_inner_nodes; i++){
    for(j = 0; j < temp_count; j++){
      if(new_nodes_flag[j] == 1){
        dist = Distance(inner_nodes_coordinate[i], temp_crackface_new_nodes_coordinate[j]);
        if(dist < layer_size){
          new_nodes_flag[i] = 0;
        }
      }
    }
  }
  for(i = 0; i < number_of_innerest_nodes; i++){
    for(j = 0; j < temp_count; j++){
      if(new_nodes_flag[j] == 1){
        dist = Distance(innerest_nodes_coordinate[i], temp_crackface_new_nodes_coordinate[j]);
        if(dist < layer_size){
          new_nodes_flag[i] = 0;
        }
      }
    }
  }
  for(i = 0; i < temp_count; i++){
    if(new_nodes_flag[i] == 1){
      AddNewNode(temp_crackface_new_nodes_coordinate[i], CRACK_FACE, OTHER);
      for(j = 0; j < temp_count; j++){
        if(new_nodes_flag[j] == 1){
          dist = Distance(temp_crackface_new_nodes_coordinate[i], temp_crackface_new_nodes_coordinate[j]);
          if(dist < layer_size){
            new_nodes_flag[j] = 0;
          }
        }
      }
    }
  }
}

void PerformCommand()
{
  ResisterCrackFrontNodes(nnodes, nodes_coordinate);
  GetPointtoPointVector(nnodes, nodes_coordinate, p_to_p_vector);
  GetUnivectoratAllPoint(nnodes, p_to_p_vector, univec_normal, univec_propa, univec_tangent);
  SetScaleFactor();
  VectorCorrectiontoSurface();
  MakeInOutLayer();

  LaminationLayer();
  LaminationPointsToNodes();
}

int main(int argc, char *argv[]){
  ReadNodes(argv[1]);//kawai.node
  ClearMinusFlag();
  ReadUnivector(argv[2]);//き裂の一番端の点の法線ベクトルの情報
  ReadPatch(argv[3]);
  sprintf(patch_file_name, "%s", argv[3]);
  ReadCrackParam(argv[4]);//param.sc_ellipse_polygon3の読み込み
  //ReadVector(argv[3]);

  PerformCommand();

  WriteCompleteNodes(argv[5]);
  WriteCrackFlag(argv[6]);
  WriteSurfaceFlag(argv[7]);
  return 0;
}
