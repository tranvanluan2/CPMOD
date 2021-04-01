/*
*      _______                       _        ____ _     _
*     |__   __|                     | |     / ____| |   | |
*        | | __ _ _ __ ___  ___  ___| |    | (___ | |___| |
*        | |/ _` | '__/ __|/ _ \/ __| |     \___ \|  ___  |
*        | | (_| | |  \__ \ (_) \__ \ |____ ____) | |   | |
*        |_|\__,_|_|  |___/\___/|___/_____/|_____/|_|   |_|
*                                                         
* -----------------------------------------------------------
*
*  TarsosLSH is developed by Joren Six at 
*  The School of Arts,
*  University College Ghent,
*  Hoogpoort 64, 9000 Ghent - Belgium
*  
* -----------------------------------------------------------
*
*  Info    : http://tarsos.0110.be/tag/TarsosLSH
*  Github  : https://github.com/JorenSix/TarsosLSH
*  Releases: http://tarsos.0110.be/releases/TarsosLSH/
* 
 */
package be.tarsos.lsh;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import mtree.tests.Data;
import mtree.utils.Constants;

/**
 * An Vector contains a vector of 'dimension' values. It serves as the main data
 * structure that is stored and retrieved. It also has an identifier (key).
 *
 * @author Joren Six
 */
public class Vector extends Data implements Serializable {

    private static final long serialVersionUID = 5169504339456492327L;

    /**
     * Values are stored here.
     */
    public Integer[] hashValues;

    public int numSucceedingNeighbors = 0;

    public int slide_range = 1;
    public double distance_to_kth_neighbor = -1;
    public int mostImportantAttribute = 0;
    public boolean checked = false;
    public boolean finishedProbeCore = false;
//    public HashMap<Integer, ArrayList<Vector>> neighborMap = new HashMap<>();
//    public ArrayList< ArrayList<ArrayList<Vector>>> neighbor_d_map = new ArrayList<>();
//    public HashMap<Integer, Integer> preceding_neighbor_map = new HashMap<>();
    public HashMap<Integer, Integer> passiveNeighbor_map = new HashMap<>();

//    public HashMap<Integer, ArrayList<Vector>> neighbor_map = new HashMap<>();
    public HashMap<Integer, Integer> neighbor_count_map = new HashMap<>();
    public HashMap<Integer, Vector> corePoints_map = new HashMap<>();
    
    public HashMap<Integer, ArrayList<Vector>> listCorePoints_map = new HashMap<>();
//    public ArrayList<Vector>[] same_slide_neighbors;
//    public ArrayList<Vector> closeNeighbor_same_slide = new ArrayList<>();
    
    public HashMap<Integer, ArrayList<Vector>> closeNeighborMap = new HashMap<>();
    public HashMap<Integer, ArrayList<Vector>> largeRNeighborCandidates = new HashMap<>();
    public HashMap<Integer, ArrayList<Vector>> linkedPointMaps = new HashMap<>();
    
    public ArrayList<Vector> neighborCores = new ArrayList<>();
    public ArrayList<Vector> closeNeighborCores = new ArrayList<>();
    public int numTotalCloseNeighbor = 0;
    public HashMap<Integer, Integer>  numCloseNeighborMap =  new HashMap<>();
//    public boolean isInCloseNeighbor = false;
    public boolean isCore = false;
    public boolean isSafe = true;
    public boolean isSmallCore = false;

//    public ArrayList<Vector> closeNeighbors 
    public ArrayList<Vector> doubleRCandidates;
//    public int num_suc_neighbors = 0;

//    public ArrayList<Vector> neighborList = new ArrayList<>();
//    public int oldestNeighborSlide = -1;
    public int lastProbRight = -1;
    public int lastProbLeft = -1;
    public int lastProbCoreLeft = -1;
    public int lastProbCoreRight = -1;
    
    public int lastProb = -1;
    public int maxChecked2RSlide = -1;
//    public ArrayList<Vector> succeedingNeighbors = new ArrayList<>();
//    public Vector closestPreNeighbor=null;
//    public double distance_to_closestPre = -1;
    public boolean added_for_probing = false;
    public boolean isInCloseNeighborList = false;

//    public HashMap<Integer, ArrayList<Vector>> numPrecedingNeighbor = new HashMap<>();
//    public int countNeighbor = 0;
    public int numNeighborPassive = 0;
    public ArrayList<Vector> closeNeighbors = new ArrayList<>();

//    public int lastProbe;
    public int sIndex;

    public double[] a;
    public double mean;
    public double std;

    public ArrayList<Integer[]> groups = new ArrayList<>();

    public double[] mean_list;
    public double[] std_list;
    public int neighborCount = 0;
    public int num_bin = 2;
    public ArrayList<Double> distance_to_neighbors = new ArrayList<>();

    public boolean checkedUB = false;
    public boolean checkedLB = false;

//    public HashMap<Integer, ArrayList<Vector>> distance_neighbor_map = new HashMap<>();
    public Vector() {

    }

    /**
     * An optional key, identifier for the vector.
     */
    private String key;

    /**
     * Creates a new vector with the requested number of dimensions.
     *
     * @param dimensions The number of dimensions.
     */
    public Vector(int dimensions) {
        this(null, new double[dimensions]);
    }
    
    public int totalNeighborCount(){
        int totalNeighbor=0;
       
        for(Integer sIdx: this.neighbor_count_map.keySet()){
            
            totalNeighbor += this.neighbor_count_map.get(sIdx);
        }
        return totalNeighbor;
    }
    
    public boolean closeToFullCore(){
        for(Vector v: closeNeighborCores){
            if(v.totalCloseNeighbor() >= Constants.k)
                return true;
        }
        return false;
    }
    
    public int totalCloseNeighbor(){
        int totalNeighbor = 0;
        for(Integer sIdx: this.linkedPointMaps.keySet()){
            totalNeighbor += this.linkedPointMaps.get(sIdx).size();
        }
        return totalNeighbor;
    }

//    public int countStoredNeighbor (){
//        int count= 0;
//        for(int sIdx: neighbor_map.keySet()){
//            count += neighbor_map.get(sIdx).size();
//        }
//        return count;
//    }
//    /**
//     * Copy constructor.
//     *
//     * @param other The other vector.
//     */
//    public Vector(Vector other) {
//        //copy the values
//        this(other.getKey(), Arrays.copyOf(other.values, other.values.length));
//    }
//    public ArrayList<Vector> getAllNeighbor() {
//        ArrayList<Vector> checked = new ArrayList<>();
//        for (ArrayList<Vector> hv : neighborMap.values()) {
//            checked.addAll(hv);
//        }
//        return checked;
//    }
    public Vector(Data d) {
        this.values = d.values;
        this.arrivalTime = d.arrivalTime;
        this.hashCode = d.hashCode;

        doubleRCandidates = new ArrayList<>();
//        distance_to_kth_neighbor = -1;

//            this.values = d.values;
//            this.hashCode = d.hashCode;
//        numPrecedingNeighbor = new HashMap<>();
        numSucceedingNeighbors = 0;
        sIndex = (int) Math.floor((arrivalTime - 1) / Constants.slide);
//        same_slide_neighbors = new ArrayList[this.num_bin];
//        for(int i =0; i < num_bin; i++){
//            this.same_slide_neighbors[i] = new ArrayList<>();
//        }
//        lastProbe = sIndex;
//        this.mean = 0;
//        for (int i = 0; i < this.values.length; i++) {
//            this.mean += this.values[i];
//        }
//        this.mean = this.mean / this.values.length;

//        this.std = 0;
//        for (int i = 0; i < this.values.length; i++) {
//            this.std += (this.values[i] - this.mean) * (this.values[i] - this.mean);
//        }
//        this.std = this.std / this.values.length;
//        this.std = Math.sqrt(this.std);
        if (null != Constants.dataFile) {
            switch (Constants.dataFile) {

                case "household2.txt":
//                    this.groups.add(new Integer[]{0, 5});
//                    this.groups.add(new Integer[]{1, 4});
//                    this.groups.add(new Integer[]{2});
//                    this.groups.add(new Integer[]{3, 6});
//                    this.groups.add(new Integer[]{0, 1, 4, 5});
//                    this.groups.add(new Integer[]{2});
//                    this.groups.add(new Integer[]{3, 6});

                    if (Constants.useLB1 || Constants.useUB1) {
                        this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6});
//this.groups.add(new Integer[]{0, 5});
//                        this.groups.add(new Integer[]{1,2,3, 4});
//                        this.groups.add(new Integer[]{2});
//                        this.groups.add(new Integer[]{2,3, 6});
//                        this.groups.add(new Integer[]{0, 1, 4, 5});
//                        this.groups.add(new Integer[]{2});
//                        this.groups.add(new Integer[]{0,5, 6});
                    }
                    break;
                case "new_hpc.txt":
//                    this.groups.add(new Integer[]{0, 5});
//                    this.groups.add(new Integer[]{1, 4});
//                    this.groups.add(new Integer[]{2});
//                    this.groups.add(new Integer[]{3, 6});
                    this.groups.add(new Integer[]{0, 1, 3});
                    this.groups.add(new Integer[]{2, 4, 5, 6});
                    break;
                case "covtype.data":
                    if (Constants.useLB1 || Constants.useUB1) {
//                        this.groups.add(new Integer[]{0});
//                        this.groups.add(new Integer[]{1});
//                        this.groups.add(new Integer[]{2});
//                        this.groups.add(new Integer[]{4});
//                        this.groups.add(new Integer[]{3});
//                        this.groups.add(new Integer[]{5});
//                        this.groups.add(new Integer[]{6});
//                        this.groups.add(new Integer[]{7});
//                        this.groups.add(new Integer[]{8});
//                        this.groups.add(new Integer[]{9});
//                        this.groups.add(new Integer[]{54});
//                        this.groups.add(new Integer[]{0,1,2,3,4,5});
//                        this.groups.add(new Integer[]{6,7,8,9,54});
                        this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                            21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                            41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54});
                    }
                    break;
                case "ethylene.txt":

                    if (Constants.useLB1 || Constants.useUB1) {

//                        this.groups.add(new Integer[]{0, 2});
////                    this.groups.add(new Integer[]{1});
//                        this.groups.add(new Integer[]{3, 4, 6, 7});
//                        this.groups.add(new Integer[]{8, 14, 15});
//                        this.groups.add(new Integer[]{1, 9});
//                        this.groups.add(new Integer[]{10, 11, 12});
//                        this.groups.add(new Integer[]{13});
                        this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13});
                    }
                    break;

                case "new_tao.txt":
                    this.groups.add(new Integer[]{0, 1, 2});
//                    this.groups.add(new Integer[]{1});
//                    this.groups.add(new Integer[]{0,2});
                    break;

                case "tao.txt":
                    if (Constants.useLB1 || Constants.useUB1) {
                        this.groups.add(new Integer[]{0, 1, 2});
                    }
//                    if (Constants.useLB2) {
//                        this.mostImportantAttribute = 1;
//                    }
//                    this.groups.add(new Integer[]{1});
//                    this.groups.add(new Integer[]{0,2});
                    break;

                case "gaussian.txt":
                    this.groups.add(new Integer[0]);
                default:
                    break;
            }
        }
        if (this.groups.size() > 0) {
            int num_group = this.groups.size();

            this.mean_list = new double[num_group];
            this.std_list = new double[num_group];
            for (int i = 0; i < num_group; i++) {
                double temp = 0;
                for (Integer j : this.groups.get(i)) {
                    temp += this.values[j];
                }
                this.mean_list[i] = temp / this.groups.get(i).length;
                temp = 0;
                for (Integer j : this.groups.get(i)) {
                    temp += (this.values[j] - this.mean_list[i]) * (this.values[j] - this.mean_list[i]);
                }

                this.std_list[i] = Math.sqrt(temp / this.groups.get(i).length);

            }
        }

    }

    /**
     * Creates a vector with the values and a key
     *
     * @param key The key of the vector.
     * @param values The values of the vector.
     */
    public Vector(String key, double[] values) {
        this.values = values;
        this.key = key;
    }

    /**
     * Moves the vector slightly, adds a value selected from -radius to +radius
     * to each element.
     *
     * @param radius The radius determines the amount to change the vector.
     */
    public void moveSlightly(double radius) {
        Random rand = new Random();
        for (int d = 0; d < getDimensions(); d++) {
            //copy the point but add or subtract a value between -radius and +radius
            double diff = radius + (-radius - radius) * rand.nextDouble();
            double point = get(d) + diff;
            set(d, point);
        }
    }

    /**
     * Set a value at a certain dimension d.
     *
     * @param dimension The dimension, index for the value.
     * @param value The value to set.
     */
    public void set(int dimension, double value) {
        values[dimension] = value;
    }

    /**
     * Returns the value at the requested dimension.
     *
     * @param dimension The dimension, index for the value.
     * @return Returns the value at the requested dimension.
     */
    public double get(int dimension) {
        return values[dimension];
    }

    /**
     * @return The number of dimensions this vector has.
     */
    public int getDimensions() {
        return values.length;
    }

    /**
     * Calculates the dot product, or scalar product, of this vector with the
     * other vector.
     *
     * @param other The other vector, should have the same number of dimensions.
     * @return The dot product of this vector with the other vector.
     * @exception ArrayIndexOutOfBoundsException when the two vectors do not
     * have the same dimensions.
     */
    public double dot(Vector other) {
        double sum = 0.0;
        for (int i = 0; i < getDimensions(); i++) {
            sum += values[i] * other.values[i];
        }
        return sum;
    }

    public void setKey(String key) {
        this.key = key;
    }

    public String getKey() {
        return key;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("values:[");
        for (int d = 0; d < getDimensions() - 1; d++) {
            sb.append(values[d]).append(",");
        }
        sb.append(values[getDimensions() - 1]).append("]");
        return sb.toString();
    }

//    public int getOldestSlide() {
//        int result = getSlideIndex();
//        for (Integer slide : numPrecedingNeighbor.keySet()) {
//            if (slide < result) {
//                result = slide;
//            }
//        }
//        return result;
//    }
//    public void addNeighbor(Vector d2) {
//        ArrayList<Vector> hv = neighborMap.get(d2.getSlideIndex());
//        if (hv == null) {
//            hv = new ArrayList<>();
//            hv.add(d2);
//            neighborMap.put(d2.getSlideIndex(), hv);
//        }
//        else{
//            hv.add(d2);
//        }
//        neighborCount += 1;
//    }
    public int getFirstSlideIdx() {
        return Math.max(getSlideIndex() - Constants.W / Constants.slide + 1, 0);
    }

//    public Integer addNeighbor_d_map(Vector d2, int bin) {
//
//        ArrayList<ArrayList<Vector>> neighbor_slide_i = this.neighbor_d_map.get(d2.getSlideIndex() - getFirstSlideIdx());
//        neighbor_slide_i.get(bin).add(d2);
//
//        return bin;
//    }
//    public Integer addNeighbor_d_map(Vector d2, double distance) {
//        Integer bin = (int) ((distance - 0.0000001) * num_bin / Constants.R);
//
//        return addNeighbor_d_map(d2, bin);
//    }
//    public void addToSameSlideNeighbor(Vector d2, int bin) {
//        this.same_slide_neighbors[bin].add(d2);
//    }
    public int countNeighbor() {

        return this.neighborCount;
    }

//    public boolean checkNeighbor(Vector d2){
//        if (DistanceFunction.LB_distance2(this, d2) <= Constants.R) {
//            
//            double d_ = DistanceFunction.UB_distance2(this, d2);
//            if(d_ <= Constants.R){
//                if ((d2.getSlideIndex() == this.getSlideIndex()-1) && (distance_to_closestPre==-1 || distance_to_closestPre > d_)){
//                    distance_to_closestPre = d_;
//                    this.closestPreNeighbor = d2;
//                }
//                return true;
//            }
//            else{
//                double exact_d = DistanceFunction.euclideanDistance(this, d2);
//                if ((d2.getSlideIndex() == this.getSlideIndex()-1) && (distance_to_closestPre==-1 || distance_to_closestPre > exact_d)){
//                    distance_to_closestPre = exact_d;
//                    this.closestPreNeighbor = d2;
//                }
//                return true;
//            }
//            
//            
//        }
//        return false;
//    }
    public int getSlideIndex() {
        return sIndex;
    }

    public boolean isOutlier() {

        return countNeighbor() < Constants.k;
    }

//        private void reset() {
//            numPrecedingNeighbor.clear();
//            numSucceedingNeighbors = 0;
//            lastProbe = sIndex;
//
//        }
//        private void updateEarliestNeighbor() {
//            
//        }
//    public boolean isUnsafeInlier() {
//
//        return !isOutlier() && numSucceedingNeighbors < Constants.k;
//    }
//    public void clean() {
//
//        numSucceedingNeighbors = 0;
//        numPrecedingNeighbor.clear();
////            clusters.clear();
//    }
//    public void clean() {
////        this.neighborList.clear();
//        this.oldestNeighborSlide =-1;
//    }
}
