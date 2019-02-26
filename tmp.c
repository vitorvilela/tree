typedef indexVector {
  int i, j, k;  
} indexVector_t;

typedef positionVector {
  double x, y, z;
} positionVector_t;

typedef box {
  int iBottom, jBottom, kBottom;
  int iTop, jTop, kTop;
} box_t;


void print_createDatatype(MPI_Comm comm, hid_t fileId, hid_t indexVectorId, hid_t positionVectorId, hid_t boxId) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  
  unsigned int attributeLenght, dataOffset, indexVectorSize, positionVectorSize, boxSize;
  int error;
  
  /*
   * Create a datatype: indexVector, positionVector and box
   */
  
  /*
   * indexVector
   */
  attributeLenght = h5tget_size(H5T_NATIVE_INTEGER, error);
  indexVectorSize = 3*attributeLenght;
  
  h5tcreate(H5T_COMPOUND, indexVectorSize, indexVectorId, error);
  
  dataOffset = 0;
  h5tinsert(iVectorId, "i", dataOffset, H5T_NATIVE_INTEGER, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(iVectorId, "j", dataOffset, H5T_NATIVE_INTEGER, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(iVectorId, "k", dataOffset, H5T_NATIVE_INTEGER, error);
  
  h5tcommit(fileId, "indexVector", indexVectorId, error);

  /*
   * positionVector
   */
  attributeLenght = h5tget_size(H5T_NATIVE_DOUBLE, error);
  positionVectorSize = 3*attributeLenght;
  
  h5tcreate(H5T_COMPOUND, positionVectorSize, positionVectorId, error);
  
  dataOffset = 0;
  h5tinsert(positionVectorId, "x", dataOffset, H5T_NATIVE_DOUBLE, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(positionVectorId, "y", dataOffset, H5T_NATIVE_DOUBLE, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(positionVectorId, "z", dataOffset, H5T_NATIVE_DOUBLE, error);
  
  h5tcommit(fileId, "positionVector", positionVectorId, error);
   
  /*
   * box
   */
  attributeLenght = h5tget_size(H5T_NATIVE_INTEGER, error);
  boxSize = 6*attributeLenght;
  
  h5tcreate(H5T_COMPOUND_F, boxSize, boxId, error);
  
  dataOffset = 0;
  h5tinsert(boxId, "iBottom", dataOffset, H5T_NATIVE_INTEGER, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(boxId, "jBottom", dataOffset, H5T_NATIVE_INTEGER, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(boxId, "kBottom", dataOffset, H5T_NATIVE_INTEGER, error); 
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(boxId, "iTop", dataOffset, H5T_NATIVE_INTEGER, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(boxId, "jTop", dataOffset, H5T_NATIVE_INTEGER, error);
  dataOffset = dataOffset + attributeLenght;
  h5tinsert(boxId, "kTop", dataOffset, H5T_NATIVE_INTEGER, error);
  
  h5tcommit(fileID, "box", boxId, error);
  
}

/*
 * Set refinment by value
 */
/*void refinment_byValue(MPI_Comm comm, leaf_t **ptr_tree, const int idxValue) {
  
  int process_rank = 0;
  MPI_Comm_rank(comm, &process_rank);
  

//   double valueMax = tree_fvMaxValue(comm, ptr_tree, idxValue);
//   double valueMin = tree_fvMinValue(comm, ptr_tree, idxValue);
//   double dValue = fabs(valueMax-valueMin);
//   double dRfn = dValue/RFN_LEVELS;


  double valueMax = 1.0;
  double valueMin = 0.0;
  double dValue = fabs(valueMax-valueMin);
  double dRfn = dValue/RFN_LEVELS;
  
  //printf("\n\nvalueMax %f, valueMin %f, dValue %f, dRfn %f\n\n", valueMax, valueMin, dValue, dRfn);
 
  double *rfnValues = malloc((RFN_LEVELS-1)*sizeof(double));
  rfnValues[0] = valueMin + dRfn; 
  //printf("\n\nrfnValues[0] %f\n", rfnValues[0]);
  for(int i=1; i<(RFN_LEVELS-1); i++) {
    rfnValues[i] = rfnValues[i-1] + dRfn; 
    //printf("rfnValues[%i] %f\n", i, rfnValues[i]);    
  }
  //printf("\n\n");
   
  for(int t=1; t<RFN_LEVELS; t++) {
  
    //printf("\nrefinment_byValue t = %i\n", t);
   
    leaf_t *iter_leaf, *tmp_leaf;
    HASH_ITER(hh, *ptr_tree, iter_leaf, tmp_leaf) {
    
      double value = iter_leaf->fv.center_vars[idxValue];
      //printf("\nvalue %f\n", value);
    
      int l = 0;     
    
      for(int i=0; i<RFN_LEVELS-1; i++) {    
	if(value > rfnValues[i]) l++;  
	else i = RFN_LEVELS;
      }
    
      if(l > (iter_leaf->current_rfn)) {      
	iter_leaf->future_rfn = iter_leaf->current_rfn + 1;
	refinment_leafRealloc(comm, ptr_tree, iter_leaf);      
      }
      else if(l < (iter_leaf->current_rfn)) {	
	iter_leaf->future_rfn = iter_leaf->current_rfn - 1;
	refinment_leafRealloc(comm, ptr_tree, iter_leaf);
      }
      
    }      
    
  }  
  
} */



void refinment_stretching(MPI_Comm comm, leaf_t **ptr_tree, const char type)

printf("\n\nvalueMax %f, valueMin %f, dValue %f, dRfn %f\n\n", valueMax, valueMin, dValue, dRfn);
printf("\n\nrfnValues[0] %f\n", rfnValues[0]);
printf("rfnValues[%i] %f\n", i, rfnValues[i]);   
printf("\n\n");
printf("\nrefinment_byValue t = %i\n", t);
printf("\nvalue %f\n", value);


void refinment_setBox(MPI_Comm comm, leaf_t **ptr_tree, double *startPoint, double *endPoint) 

printf("\nif(leaf_xLB %f >= startPoint[0] %f && leaf_yLB %f >= startPoint[1] %f && leaf_zLB %f >= startPoint[2] %f && leaf_xUB %f <= endPoint[0] %f && leaf_yUB %f <= endPoint[1] %f && leaf_zUB %f <= endPoint[2] %f\n", leaf_xLB, startPoint[0], leaf_yLB, startPoint[1], leaf_zLB, startPoint[2], leaf_xUB, endPoint[0], leaf_yUB, endPoint[1], leaf_zUB, endPoint[2]);


  fprintf(ptr_monitor, "\n\t\t \n");
  
  
      //DBG
    fclose(*ptr_monitor);