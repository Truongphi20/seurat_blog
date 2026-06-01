


calculate_mu <- function(Beta, model_matrix, offset_matrix){
  make_offset_hdf5_mat <- is(offset_matrix, "DelayedMatrix")
  if(make_offset_hdf5_mat){
    mu <- exp(delayed_matrix_multiply(DelayedArray::DelayedArray(Beta),
                                      DelayedArray::DelayedArray(t(model_matrix))) + offset_matrix)
    mu <- HDF5Array::writeHDF5Array(mu)
    mu
  }else{
    exp(Beta %*% t(model_matrix) + offset_matrix)
  }
}
