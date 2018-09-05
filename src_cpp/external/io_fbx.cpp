#include <vector>
#include <string>
#include <iostream>
#include <map>

#include "fbxsdk.h"

#include "delfem2/../../external/io_fbx.h"

class CVal3
{
public:
  CVal3(double v0,double v1,double v2){
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
  }
public:
  double v0,v1,v2;
};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Recursive_Mesh( FbxNode *pNode, std::vector<FbxMesh*>& apMesh)
{
  FbxNodeAttribute *pAttrib = pNode->GetNodeAttribute();
  if( pAttrib ){
    FbxNodeAttribute::EType type = pAttrib->GetAttributeType();
    if( type == FbxNodeAttribute::eMesh ){
      FbxMesh *pMesh = (FbxMesh*)pAttrib;
      apMesh.push_back(pMesh);
    }
  }
  //////
  int childNodeNum = pNode->GetChildCount();
  for ( int ichild = 0; ichild < childNodeNum; ichild++ )
  {
    FbxNode *pChild = pNode->GetChild(ichild);  // 子ノードを取得
    Recursive_Mesh( pChild,apMesh);
  }
}

void Recursive_Skeleton
( FbxNode *pNode,
 std::vector< std::pair<std::string,std::string> >& aSkeletonName)
{
  FbxNodeAttribute *pAttrib = pNode->GetNodeAttribute();
  if( pAttrib ){
    FbxNodeAttribute::EType type = pAttrib->GetAttributeType();
    if( type == FbxNodeAttribute::eSkeleton ){
      ////
      std::string name = pNode->GetName();
      std::string name_p = "";
      {
        FbxNodeAttribute *pAttrib_p = pNode->GetParent()->GetNodeAttribute();
        FbxNodeAttribute::EType type_p = FbxNodeAttribute::eNull;
        if( pAttrib_p ){ type_p = pAttrib_p->GetAttributeType(); }
        if( type_p  == FbxNodeAttribute::eSkeleton ){
          name_p = pNode->GetParent()->GetName();
        }
      }
      aSkeletonName.push_back( std::make_pair(name,name_p) );
    }
  }
  //////
  int childNodeNum = pNode->GetChildCount();
  for ( int ichild = 0; ichild < childNodeNum; ichild++ ){
    FbxNode *pChild = pNode->GetChild(ichild);  // 子ノードを取得
    Recursive_Skeleton(pChild,aSkeletonName);
  }
}

void getXYZElem_Mesh
(const FbxMesh* pMesh,
 std::vector<double>& aXYZ,
 std::vector<int>& aElemInd,
 std::vector<int>& aElem)
{
  aElem.clear();
  aElemInd.clear();
  {
    const int nelemind = pMesh->GetPolygonVertexCount();
    {
      const int* pind = pMesh->GetPolygonVertices();
      aElem.assign(pind,pind+nelemind);
    }
    const int nelem = pMesh->GetPolygonCount();
    aElemInd.assign(nelem+1,0);
    for ( int ielem = 0; ielem < nelem; ielem++ ) {
      int npoel = pMesh->GetPolygonSize( ielem );
      aElemInd[ielem+1] += npoel;
    }
    for( int ielem = 0;ielem<nelem;++ielem){
      aElemInd[ielem+1] += aElemInd[ielem];
    }
    assert( nelemind == aElemInd[nelem] );
  }
  ////
  double offset[3] = {0.0, 0.0, 0.0};
  double scale=1;
  {
    const FbxNode* pNode = pMesh->GetNode();
    FbxVector4 vec4 = pNode->LclTranslation.Get();
    FbxVectorTemplate3<double> scale3 = pNode->LclScaling.Get();
    {
      offset[0] = vec4[0];
      offset[1] = vec4[1];
      offset[2] = vec4[2];
    }
//    std::cout << vec4[0] << " " << vec4[1] << " " << vec4[2] << std::endl;
//    std::cout << scale3[0] << " " << scale3[1] << " " << scale3[2] << std::endl;
//    scale = scale3[0];
  }
//  std::cout << "scale: " << scale << std::endl;
  aXYZ.clear();
  {
    const int npoint = pMesh->GetControlPointsCount();
    const FbxVector4* src = pMesh->GetControlPoints();
    for ( int ipoint = 0; ipoint < npoint; ++ipoint ) {
      aXYZ.push_back(src[ipoint][0]*scale+offset[0]);
      aXYZ.push_back(src[ipoint][1]*scale+offset[1]);
      aXYZ.push_back(src[ipoint][2]*scale+offset[2]);
    }
  }
}

void GetMaterialProperty
(const FbxSurfaceMaterial * pMaterial,
 const char * pPropertyName,
 const char * pFactorPropertyName,
 double aRGB[3],
 std::vector<CTextureInfo_RigMsh>& aMatTex)
{
  const FbxProperty material_propaty = pMaterial->FindProperty(pPropertyName);
  if (!material_propaty.IsValid()){ return; }
  {
    FbxDouble3 lResult = material_propaty.Get<FbxDouble3>();
    aRGB[0] = lResult[0];
    aRGB[1] = lResult[1];
    aRGB[2] = lResult[2];
    const FbxProperty lFactorProperty = pMaterial->FindProperty(pFactorPropertyName);
    if (lFactorProperty.IsValid()){
      double lFactor = lFactorProperty.Get<FbxDouble>();
      aRGB[0] *= lFactor;
      aRGB[1] *= lFactor;
      aRGB[2] *= lFactor;
    }
  }
  /*
  {
    const FbxProperty mpc = material_propaty.GetDstProperty();
    std::cout << mpc.GetName() << std::endl;
    const int ntexture = mpc.GetSrcObjectCount<FbxFileTexture>();
    for(int itexture=0;itexture<ntexture;itexture++){
      const FbxFileTexture* pFileTex = mpc.GetSrcObject<FbxFileTexture>(itexture);
      if(!pFileTex){ continue; }
      std::cout << "parent tex " << pFileTex->GetFileName() << std::endl;
    }
  }
   */
  {
    aMatTex.clear();
    const int ntexture = material_propaty.GetSrcObjectCount<FbxFileTexture>();
    for(int itexture=0;itexture<ntexture;itexture++){
      const FbxFileTexture* pFileTex = material_propaty.GetSrcObject<FbxFileTexture>(itexture);
      if(!pFileTex){ continue; }
      CTextureInfo_RigMsh mat_tex;
//      mat_tex.full_path = pFileTex->GetFileName();
      mat_tex.full_path = pFileTex->GetFileName();
      mat_tex.name = pFileTex->GetName();
      mat_tex.uv_setname = pFileTex->UVSet.Get();
      aMatTex.push_back(mat_tex);
    }
    if( aMatTex.size() == 0 ){
      std::cout << "no texture this material" << std::endl;
    }
  }
  
/*
  {
    const int nLayeredTex = material_propaty.GetSrcObjectCount<FbxLayeredTexture>();
    for(int iLayeredTex=0;iLayeredTex<nLayeredTex;iLayeredTex++){
      FbxLayeredTexture* lLayeredTexture = material_propaty.GetSrcObject<FbxLayeredTexture>(iLayeredTex);
      std::cout << "  ilayerdTexture" << iLayeredTex << std::endl;
      const int nTexFile = lLayeredTexture->GetSrcObjectCount<FbxFileTexture>();
      for(int j=0;j<nTexFile;j++){
        FbxFileTexture* lFileTexture = lLayeredTexture->GetSrcObject<FbxFileTexture>(j);
        if(!lFileTexture){ continue; }
        FbxString uvsetName = lFileTexture->UVSet.Get();
        std::string uvSetString = uvsetName.Buffer();
        std::string filepath = lFileTexture->GetFileName();
        std::cout << "  hugeahgea:" << uvsetName << " " << uvSetString << " " << filepath << std::endl;
      }
    }
  }
  */
}

void getMaterial(FbxMesh* pMesh,std::vector<CMaterial_RigMsh>& aMaterial)
{
  FbxNode* node = pMesh->GetNode();
  if ( node == 0 ) { return; }
  const int nmaterial = node->GetMaterialCount();
  aMaterial.resize(nmaterial);
  for( int imaterial = 0; imaterial < nmaterial; ++imaterial ) {
    const FbxSurfaceMaterial* pSurfMat = node->GetMaterial( imaterial );
    if ( pSurfMat == 0 ) { continue; }
    GetMaterialProperty(pSurfMat,
                        FbxSurfaceMaterial::sDiffuse,
                        FbxSurfaceMaterial::sDiffuseFactor,
                        aMaterial[imaterial].RGB_Diffuse,
                        aMaterial[imaterial].aTexture_Diffuse);
  }
}



void getMaterialMapping_Layer
(const FbxLayer* pLayer,
 std::string& material_mapping_mode,
 std::vector<int>& material_mapping_face)
{
  const FbxLayerElementMaterial* pLayerElemMat = pLayer->GetMaterials();
  if (!pLayerElemMat){ return; }
  FbxLayerElement::EMappingMode mappingMode = pLayerElemMat->GetMappingMode();
  if( mappingMode == FbxLayerElement::eAllSame ){
    material_mapping_mode = "ALL_SAME";
  }
  else if( mappingMode == FbxLayerElement::eByPolygon ){
    material_mapping_mode = "BY_FACE";
    if (pLayerElemMat->GetReferenceMode() == FbxLayerElement::eIndexToDirect ){
      const int nelem = pLayerElemMat->GetIndexArray().GetCount();
      material_mapping_face.resize(nelem);
      for (int ielem = 0; ielem < nelem; ielem++) {
        int imaterial = pLayerElemMat->GetIndexArray().GetAt(ielem);
        material_mapping_face[ielem] = imaterial;
      }
    }
    else{
      std::cout << "unsupported ref mode" << std::endl;
    }
  }
  else{
    std::cout << "unknown mapping mode: " << mappingMode << std::endl;
  }
}

/*
 void LoadUVInformation2(FbxMesh* pMesh)
 {
 FbxStringList  uvsetName;
 pMesh->GetUVSetNames(uvsetName);
 int numUVSet = uvsetName.GetCount();
 std::cout << numUVSet << std::endl;
 
 std::cout << "GetElementUVCount: " << pMesh->GetElementUVCount() << std::endl;
 
 bool unmapped = false;
 //  FbxGeometryElementUV* leUV = pMesh->GetElementUV(l);
 
 int lPolygonCount = pMesh->GetPolygonCount();
 for(int uv=0;uv<numUVSet;uv++)
 {
 //    meshNode->uvsetID[uvsetName.GetStringAt(uv)] = uv;
 for(int ipoly=0;ipoly<lPolygonCount;ipoly++)
 {
 int lPolygonsize = pMesh->GetPolygonSize(ipoly);
 for(int pol=0;pol<lPolygonsize;pol++)
 {
 int itxuvind = pMesh->GetTextureUVIndex(ipoly,pol);
 //        leUV->GetDirectArray().GetAt(itxuvind);
 FbxString name = uvsetName.GetStringAt(uv);
 FbxVector2 texCoord;
 pMesh->GetPolygonVertexUV( ipoly, pol, name, texCoord, unmapped);
 //        std::cout << ipoly << " " << pol << " " << name << " " << itxuvind << " " << texCoord[0] << " " << texCoord[1] << " " << unmapped << std::endl;
 //        meshNode->m_texcoordArray.push_back(texCoord);
 }
 }
 }
 }
 */

void getUV_Geom(FbxMesh* pMesh)
{
  //get all UV set names
  FbxStringList listUVSetName;
  pMesh->GetUVSetNames(listUVSetName);
  std::cout << "number of uv setname" << listUVSetName.GetCount() << std::endl;
  for (int iUVSet = 0; iUVSet < listUVSetName.GetCount(); iUVSet++){ //iterating over all uv sets
    const char* sUVSetName = listUVSetName.GetStringAt(iUVSet);
    std::cout << "  uv  setname: " << sUVSetName << std::endl;
    const FbxGeometryElementUV* pGeomElemUV = pMesh->GetElementUV(sUVSetName);
    if(!pGeomElemUV){ continue; }
    
    // only support mapping mode eByPolygonVertex and eByControlPoint
    if( pGeomElemUV->GetMappingMode() != FbxGeometryElement::eByPolygonVertex &&
       pGeomElemUV->GetMappingMode() != FbxGeometryElement::eByControlPoint )
      return;
    
    const bool lUseIndex = pGeomElemUV->GetReferenceMode() != FbxGeometryElement::eDirect;
    const int lIndexCount= (lUseIndex) ? pGeomElemUV->GetIndexArray().GetCount() : 0;
    
    const int nelem = pMesh->GetPolygonCount();
//    std::cout << "  nface: " << nelem << " " << "npoint: " << lIndexCount << std::endl;
    
    if( pGeomElemUV->GetMappingMode() == FbxGeometryElement::eByControlPoint ){
      for( int ielem = 0; ielem < nelem; ++ielem ){
        const int npoel = pMesh->GetPolygonSize(ielem);
        for( int ipoel = 0; ipoel < npoel; ++ipoel ){
          const int ipoint = pMesh->GetPolygonVertex(ielem,ipoel);
          const int iUV = lUseIndex ? pGeomElemUV->GetIndexArray().GetAt(ipoint) : ipoint;
          FbxVector2 lUVValue = pGeomElemUV->GetDirectArray().GetAt(iUV);
        }
      }
    }
    else if (pGeomElemUV->GetMappingMode() == FbxGeometryElement::eByPolygonVertex){
      int lPolyIndexCounter = 0;
      for( int ielem = 0; ielem < nelem; ++ielem ){
        const int npoel = pMesh->GetPolygonSize(ielem);
        for( int ipoel = 0; ipoel < npoel; ++ipoel ){
          if (lPolyIndexCounter >= lIndexCount){ break; }
          //the UV index depends on the reference mode
          const int iUV = lUseIndex ? pGeomElemUV->GetIndexArray().GetAt(lPolyIndexCounter) : lPolyIndexCounter;
          FbxVector2 lUVValue = pGeomElemUV->GetDirectArray().GetAt(iUV);
          lPolyIndexCounter++;
        }
      }
    }
    ///////
    for(int iface=0;iface<nelem;iface++){
      int npofa = pMesh->GetPolygonSize(iface);
      for(int ipofa=0;ipofa<npofa;ipofa++){
        //        int itxuvind = pMesh->GetTextureUVIndex(iface,ipofa);
        FbxVector2 texCoord;
        bool unmapped;
        pMesh->GetPolygonVertexUV( iface, ipofa, sUVSetName, texCoord, unmapped);
        //        std::cout << ipoly << " " << pol << " " << name << " " << itxuvind << " " << texCoord[0] << " " << texCoord[1] << " " << unmapped << std::endl;
      }
    }
  }
}


void getUV_Layer
(const FbxLayer* pLayer,
 std::string& uv_setname,
 std::string& uv_bindmode,
 std::vector<double>& aUV)
{
  const FbxLayerElementUV* pLayerElemUV = pLayer->GetUVs();
  if ( pLayerElemUV == 0 ) { return; }
  uv_setname = pLayerElemUV->GetName();
  //////
  FbxLayerElement::EMappingMode mappingMode = pLayerElemUV->GetMappingMode();
  if ( mappingMode == FbxLayerElement::eByPolygonVertex ) {
    uv_bindmode = "POEL";
    FbxLayerElement::EReferenceMode refMode     = pLayerElemUV->GetReferenceMode();
    if ( refMode == FbxLayerElement::eDirect ) {
      const int nUV  = pLayerElemUV->GetDirectArray().GetCount();
      aUV.resize(nUV*2);
      for ( int iUV = 0; iUV < nUV; ++iUV ) {
        aUV[iUV*2+0] = (double)pLayerElemUV->GetDirectArray().GetAt( iUV )[ 0 ];
        aUV[iUV*2+1] = (double)pLayerElemUV->GetDirectArray().GetAt( iUV )[ 1 ];
      }
    }
    else if ( refMode == FbxLayerElement::eIndexToDirect ) {
      const int nUV  = pLayerElemUV->GetIndexArray().GetCount();
      aUV.resize(nUV*2);
      for ( int iUV = 0; iUV < nUV; ++iUV ) {
        int index = pLayerElemUV->GetIndexArray().GetAt( iUV );
        aUV[iUV*2+0] = (double)pLayerElemUV->GetDirectArray().GetAt( index )[ 0 ];
        aUV[iUV*2+1] = (double)pLayerElemUV->GetDirectArray().GetAt( index )[ 1 ];
      }
    }
  }
  else if ( mappingMode == FbxLayerElement::eByControlPoint){
  }
}


void getRigging
(const FbxMesh* pMesh,
 std::vector<CSkin_RigMsh>& aSkin,
 std::vector<std::vector<CVal3> >& aPosSkel,
 const std::map<std::string,int>& mapName2Indb)
{
  const int nskin  = pMesh->GetDeformerCount( FbxDeformer::eSkin );
  double offset_mesh[3] = {0.0, 0.0, 0.0};
  {
    const FbxNode* pNode = pMesh->GetNode();
    FbxVector4 vec4 = pNode->LclTranslation.Get();
    offset_mesh[0] = vec4[0];
    offset_mesh[1] = vec4[1];
    offset_mesh[2] = vec4[2];
  }
  aSkin.resize(nskin);
  for ( int iskin = 0; iskin < nskin; ++iskin ) {
    const FbxSkin* pSkin = (FbxSkin*)(pMesh->GetDeformer( iskin, FbxDeformer::eSkin ));
    const int nbone = pSkin->GetClusterCount();
    aSkin[iskin].aBone.resize(nbone);
    for (int ibone = 0; ibone < nbone; ++ibone ) {
      const FbxCluster* pCluster = pSkin->GetCluster( ibone );
      {
        const FbxNode* pNode = pCluster->GetLink();
        if( pNode != 0 ){
          aSkin[iskin].aBone[ibone].name = pNode->GetName();
          const FbxNodeAttribute *pAttrib = pNode->GetNodeAttribute();
          assert( pAttrib != 0 );
          FbxNodeAttribute::EType type = pAttrib->GetAttributeType();
          assert( type == FbxNodeAttribute::eSkeleton );
        }
      }
      { // get initial pos
        FbxAMatrix matLink;  pCluster->GetTransformLinkMatrix( matLink  );
        FbxAMatrix matTrans; pCluster->GetTransformMatrix(     matTrans );
        FbxVector4 posLink = matLink.GetT();
        double pos[3] = {posLink[0], posLink[1], posLink[2]};
        double pos1[3] = {0,0,0};
        for(int i=0;i<3;++i){
          for(int j=0;j<3;++j){
            pos1[i] += matTrans.Get(i, j)*pos[j];
          }
        }
        {
          FbxVector4 scale = matTrans.GetS();
//          std::cout << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
        }
//        pos1[0] += matTrans.GetT()[0];
//        pos1[1] += matTrans.GetT()[1];
//        pos1[2] += matTrans.GetT()[2];
        pos1[0] += offset_mesh[0];
        pos1[1] += offset_mesh[1];
        pos1[2] += offset_mesh[2];
//        std::cout << pos1[0] << " " << pos1[1] << " " << pos1[2] << std::endl;
        const std::string name = aSkin[iskin].aBone[ibone].name;
        const std::map<std::string,int>::const_iterator itr = mapName2Indb.find(name);
        if( itr != mapName2Indb.end() ){
          int iskel = itr->second;
          CVal3 v(pos1[0],pos1[1],pos1[2]);
          aPosSkel[iskel].push_back(v);
        }
      }
      { // get index weight
        const int  npocl = pCluster->GetControlPointIndicesCount();
        aSkin[iskin].aBone[ibone].aIpWeight.resize(npocl);
        const int*  pointAry = pCluster->GetControlPointIndices();
        const double* weightAry = pCluster->GetControlPointWeights();
        for ( int ipocl = 0; ipocl < npocl; ++ipocl ) {
          const int ipoint  = pointAry[ ipocl ];
          const double weight = (double)weightAry[ ipocl ];
          aSkin[iskin].aBone[ibone].aIpWeight[ipocl] = std::make_pair(ipoint,weight);
        }
      }
    }
  }
}

void Read_FBX(const std::string& path, CRigMsh& rigmsh)
{
  ////
  FbxNode *pRootNode = 0;
  {
    FbxManager *lSdkManager = FbxManager::Create();
    FbxIOSettings * ios = FbxIOSettings::Create(lSdkManager, IOSROOT );// create an IOSettings object
    // set some IOSettings options
    ios->SetBoolProp(IMP_FBX_MATERIAL, true);
    ios->SetBoolProp(IMP_FBX_TEXTURE,  true);
    FbxScene* lScene = FbxScene::Create(lSdkManager,""); // create an empty scene
    {
      FbxImporter *lImporter = FbxImporter::Create( lSdkManager, "");
      lImporter->Initialize(path.c_str(), -1, ios);
      lImporter->Import(lScene); // Import the scene.
      lImporter->Destroy(); // Destroy the importer.
    }
    pRootNode = lScene->GetRootNode();
  }
  ////////////////////////////////
  rigmsh.Clear();
  std::vector<FbxMesh*> apMesh;
  { // initialize mesh
    Recursive_Mesh(pRootNode,apMesh);
    rigmsh.aMesh.resize( apMesh.size() );
  }
  std::map<std::string,int> mapName2Indb;
  { // initialize skeleton
    std::vector< std::pair<std::string,std::string> > aBoneName;
    Recursive_Skeleton(pRootNode,aBoneName);
    const int nbone = aBoneName.size();
    rigmsh.aBone.resize(nbone);
    for(int ibone=0;ibone<nbone;++ibone){
      rigmsh.aBone[ibone].name = aBoneName[ibone].first;
    }
    mapName2Indb.clear();
    for(int ibone=0;ibone<nbone;++ibone){
      mapName2Indb.insert( std::make_pair(rigmsh.aBone[ibone].name,ibone) );
    }
    for(int ibone=0;ibone<nbone;++ibone){
      std::map<std::string,int>::iterator itr_p = mapName2Indb.find(aBoneName[ibone].second);
      int iskel_p = (itr_p==mapName2Indb.end()) ? -1:itr_p->second;
      rigmsh.aBone[ibone].ibone_parent = iskel_p;
    }
  }
  std::vector< std::vector<CVal3> > aPosBone(rigmsh.aBone.size());
  for(int imesh=0;imesh<apMesh.size();++imesh){
    CMesh_RigMsh& mesh = rigmsh.aMesh[imesh];
    getXYZElem_Mesh(apMesh[imesh],mesh.aXYZ_ini,mesh.aElemInd,mesh.aElem);
    mesh.aXYZ = mesh.aXYZ_ini;
    //////////////
    const int nlayer = apMesh[imesh]->GetLayerCount();
    rigmsh.aMesh[imesh].aLayer.resize(nlayer);
    for ( int ilayer = 0; ilayer < nlayer; ++ilayer ) {
      const FbxLayer* pLayer = apMesh[imesh]->GetLayer( ilayer );
      CLayer_RigMsh& layer = mesh.aLayer[ilayer];
      getUV_Layer(pLayer,
                  layer.uv_setname,
                  layer.uv_bindmode,
                  layer.aUV);
      getMaterialMapping_Layer(pLayer,
                               layer.material_mapping_mode,
                               layer.material_mapping_face);
      layer.Finalize();
    }
    /////
    getMaterial(apMesh[imesh],mesh.aMaterial);
    { // rigging information
      getRigging(apMesh[imesh],
                 mesh.aSkin, aPosBone, mapName2Indb);
      for(int iskin=0;iskin<mesh.aSkin.size();++iskin){
        mesh.aSkin[iskin].Finalize(mesh.aXYZ_ini.size()/3);
      }
    }
  }
  for(int ibone=0;ibone<aPosBone.size();++ibone){
    double x=0.0, y=0.0, z=0.0;
    for(int iv=0;iv<aPosBone[ibone].size();++iv){
      x += aPosBone[ibone][iv].v0;
      y += aPosBone[ibone][iv].v1;
      z += aPosBone[ibone][iv].v2;
//      std::cout << ibone << "   (" << aPosBone[ibone][iv].v0 << " " << aPosBone[ibone][iv].v1 << " " << aPosBone[ibone][iv].v2 << ") " << std::endl;
    }
    if( aPosBone[ibone].size() > 0 ){
      rigmsh.aBone[ibone].is_active = true;
//      rigmsh.aBone[ibone].pos_ini[0] = x/aPosBone[ibone].size();
//      rigmsh.aBone[ibone].pos_ini[1] = y/aPosBone[ibone].size();
//      rigmsh.aBone[ibone].pos_ini[2] = z/aPosBone[ibone].size();
      rigmsh.aBone[ibone].pos_ini[0] = aPosBone[ibone][0].v0;
      rigmsh.aBone[ibone].pos_ini[1] = aPosBone[ibone][0].v1;
      rigmsh.aBone[ibone].pos_ini[2] = aPosBone[ibone][0].v2;
    }
    else{
      rigmsh.aBone[ibone].is_active = false;
    }
    rigmsh.aBone[ibone].pos[0] = rigmsh.aBone[ibone].pos_ini[0];
    rigmsh.aBone[ibone].pos[1] = rigmsh.aBone[ibone].pos_ini[1];
    rigmsh.aBone[ibone].pos[2] = rigmsh.aBone[ibone].pos_ini[2];
    rigmsh.aBone[ibone].quat_joint[0] = 1;
    rigmsh.aBone[ibone].quat_joint[1] = 0;
    rigmsh.aBone[ibone].quat_joint[2] = 0;
    rigmsh.aBone[ibone].quat_joint[3] = 0;
  }
  rigmsh.UpdateBonePos();
}
