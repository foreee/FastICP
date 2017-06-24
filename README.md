# FastICP (point to point)
 ICP algorithm , fast and easy to use in point cloud registration . 
 
 Details are found in:
 
 Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes. 
 IEEE Transactions on pattern analysis and machine intelligence, 239-256.
 
 The speed of demo under release X86 mode is close to the funtion "pcregrigid" in Matlab2015a.
 
 Nearest points searching is based on libnabo.from http://github.com/ethz-asl/libnabo.
 
 Include the Extrapolation in quaternion space. Details are found in:
 
 Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes. 
 IEEE Transactions on pattern analysis and machine intelligence, 239-256.
 
 # Tutorial
 
 1- Download Eigen 3.0+ from http://eigen.tuxfamily.org.
 
    Add the "Eigen" to fold "src/nabo/".
    
 2- Download boost form http://www.boost.org.
 
    Set path of boost to your project.
    
 3- Use FastICP as the file "src/main.c" shown.
 
 
 



