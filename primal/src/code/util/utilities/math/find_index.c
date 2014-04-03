long find_index(long val, long *data, long length){
    long ans, i;
    for(i=0;i<length;i++){
        if (data[i] == val)
            return(i);
    }
    return(-999);
}
