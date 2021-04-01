/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import mtree.tests.Data;
import mtree.tests.MTTest;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class ApproxCPOD2 {

    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public static HashMap<Integer, ArrayList<CorePoint>> all_core_points = new HashMap<>();
//    public static HashMap<Integer, MTreeCorePoint> all_indexed_cores = new HashMap<>();
    public static MTreeCorePoint mtree = new MTreeCorePoint();

    public static ArrayList<CorePoint> all_distinct_cores = new ArrayList<>();
    public static HashMap<Integer, HashSet<C_Data>> outlierList = new HashMap<>();

//    public static HashMap<Integer, HashSet<CorePoint>> corePointMapTrigger = new HashMap<>();
    public static HashMap<Integer, HashSet<C_Data>> neighborCountTrigger = new HashMap<>();
//    public static HashMap<Integer, 
//    public double timeFindingCore = 0;
//    public double timeCheckingCandidates = 0;
//    public double timeAddingToNeighborCount = 0;
//    public double timeMarkProbing = 0;
//    public double avg_candidate_size = 0;
//    public double timeForAddingToOutlierList = 0;

    public static double timeProcessExpiredSlide = 0;
    public static double timeCreatingCore = 0;
    public static double timeProbing = 0;
    public static double timeReProbing = 0;
//    public static int count = 0;
    public static double numPointNeedProb = 0;
    public static double avg_points_check = 0;
    public static int countPoint = 0;
//    public static long numDCS = 0;
//    public static long numPointNeedNS = 0;
//    public static long numDCS = 0;
//    public static long count = 0;
//    public static long numDCSForIndexing = 0;

    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {

        currentTime = _currentTime;

        ArrayList<C_Data> d_to_process = new ArrayList<>(data.size());
//        HashSet<Integer> slide_to_process = new HashSet<>();
//        System.out.println("Current time = "+ currentTime);
        int[] slide_to_process;
        if (currentTime == W) {
            slide_to_process = new int[W / slide];
            for (int i = 0; i < slide_to_process.length; i++) {
                slide_to_process[i] = i;
            }
        } else {
            MTTest.count += 1;
            slide_to_process = new int[]{(currentTime - 1) / slide};
        }

        for (int i = 0; i < data.size(); i++) {
            Data o = data.get(i);
            C_Data d = new C_Data(o);

            d_to_process.add(d);

            ArrayList<C_Data> idx = all_slides.get(d.sIndex);
            if (idx != null) {

                idx.add(d);
            } else {
                idx = new ArrayList<>(Constants.slide);
                idx.add(d);
                all_slides.put(d.sIndex, idx);
            }
        }
        MTTest.start = Utils.getCPUTime();

//        long start = Utils.getCPUTime();
        ArrayList<Data> result = new ArrayList<>((int) (data.size() * 1.0 / 100));
        expiredSlideIndex = (currentTime - 1) / slide - Constants.W / slide;
//        System.out.println("Expire slide index = " + expiredSlideIndex);
//        long start = Utils.getCPUTime();
        processExpiredData(expiredSlideIndex);
//        System.out.println("Selecting core!");

        for (int sIdx : slide_to_process) {

            ArrayList<CorePoint> corePoints = selectCore(sIdx);

            all_core_points.put(sIdx, corePoints);
//            all_indexed_cores.put(sIdx, mtree);
//            System.out.println("Num new core = " + corePoints.size());
//            System.out.println("All distinc core = " + all_distinct_cores.size());
        }
//        System.out.println("Finished selecting cores!");

        int newestSlide = (currentTime - 1) / Constants.slide;

        if (currentTime == Constants.W) {
            for (CorePoint c : all_distinct_cores) {
                c.totalHalfRPoints = c.getTotalHalfRPoints();
            }
        } else if (data.size() == Constants.slide) {
            for (CorePoint c : all_distinct_cores) {
                if (c.closeNeighbors_halfR.get(newestSlide) != null) {
                    c.totalHalfRPoints += c.closeNeighbors_halfR.get(newestSlide).size();
                }
            }
        }
//        System.out.println("Probing current slide");
        for (int i = 0; i < d_to_process.size(); i++) {
            C_Data d = d_to_process.get(i);
            d.neighborCount = probe_current_slide(d);
            d.numSucceedingNeighbor = d.neighborCount;
        }
//        System.out.println("Probing!");

        for (Integer sidx : all_slides.keySet()) {
            for (int i = 0; i < all_slides.get(sidx).size(); i++) {
                C_Data d = all_slides.get(sidx).get(i);
                if (d.neighborCount < Constants.k) {
//                    System.out.println("Probing for d.arrivaltime=" + d.arrivalTime);
                    prob(d, newestSlide);
                }
                if (d.neighborCount < Constants.k) {
                    result.add(d);
                }
            }
        }

        return result;

    }

    private void prob(C_Data d, int newestSlide) {
        if (d.numSucceedingNeighbor >= Constants.k) {
            return;
        }
        if (d.closeCoreMaps_halfR != null && d.closeCoreMaps_halfR.totalHalfRPoints >= Constants.k + 1) {
            return;
        }

        //probing right
        for (int sidx = d.lastProbRight + 1; sidx <= newestSlide; sidx++) {
            prob_slide(d, sidx);
            if (d.neighborCount >= Constants.k) {
                break;
            }

        }
        //prob left
        if (d.neighborCount < Constants.k) {
            for (int sidx = d.lastProbLeft - 1; sidx > newestSlide - Constants.W / Constants.slide; sidx--) {
                prob_slide(d, sidx);
                if (d.neighborCount >= Constants.k) {
                    break;
                }
            }
        }
    }

    private int prob_slide(C_Data d, int sidx) {
        int numNeighbors;
        ArrayList<CorePoint> corePoints = all_core_points.get(sidx);
        if (d.closeCoreMaps_halfR != null && corePoints.contains(d.closeCoreMaps_halfR)) {
            int candidate_size = d.closeCoreMaps_halfR.closeNeighbors_halfR.get(sidx).size()
                    + d.closeCoreMaps_halfR.closeNeighbors_R.get(sidx).size()
                    + d.closeCoreMaps_halfR.closeNeighbors_3halfR.get(sidx).size();
            int candidate_same_slide = d.closeCoreMaps_halfR.closeNeighbors_halfR.get(d.sIndex).size()
                    + d.closeCoreMaps_halfR.closeNeighbors_R.get(d.sIndex).size()
                    + d.closeCoreMaps_halfR.closeNeighbors_3halfR.get(d.sIndex).size();
            numNeighbors = d.numNeighborSameSlide * candidate_size / candidate_same_slide;

        } else {
            ArrayList<CorePoint> commonCorePoints = new ArrayList<>();
            for (CorePoint c : d.closeCoreMaps_R) {
                if (corePoints.contains(c)) {
                    commonCorePoints.add(c);
                }
            }
//            System.out.println("Num #common cores= " + commonCorePoints.size());
            if (!commonCorePoints.isEmpty()) {
                CorePoint c = commonCorePoints.get(0);
                int candidate_size = 0;

                candidate_size += c.closeNeighbors_halfR.get(sidx).size() + c.closeNeighbors_R.get(sidx).size()
                        + c.closeNeighbors_3halfR.get(sidx).size() + c.closeNeighbors_2R.get(sidx).size();

                int candidate_same_slide = 0;

                candidate_same_slide += c.closeNeighbors_halfR.get(d.sIndex).size() + c.closeNeighbors_R.get(d.sIndex).size()
                        + c.closeNeighbors_3halfR.get(d.sIndex).size() + c.closeNeighbors_2R.get(d.sIndex).size();

                numNeighbors = d.numNeighborSameSlide * candidate_size / candidate_same_slide;
            } else {
                for (CorePoint c : d.closeCoreMaps_2R) {
                    if (corePoints.contains(c)) {
                        commonCorePoints.add(c);
                    }
                }

                if (!commonCorePoints.isEmpty()) {
                    int candidate_size = 0;
                    for (CorePoint c : commonCorePoints) {
                        candidate_size += c.closeNeighbors_halfR.get(sidx).size() + c.closeNeighbors_R.get(sidx).size()
                                + c.closeNeighbors_3halfR.get(sidx).size() + c.closeNeighbors_2R.get(sidx).size();

                    }
                    int candidate_same_slide = 0;
                    for (CorePoint c : commonCorePoints) {
                        candidate_same_slide += c.closeNeighbors_halfR.get(d.sIndex).size() + c.closeNeighbors_R.get(d.sIndex).size()
                                + c.closeNeighbors_3halfR.get(d.sIndex).size() + c.closeNeighbors_2R.get(d.sIndex).size();

                    }

                    numNeighbors = d.numNeighborSameSlide * candidate_size / candidate_same_slide;
                } else {
                    numNeighbors = 0;
                }

            }
        }

        d.neighborCount += numNeighbors;
        if (sidx < d.sIndex) {
            d.lastProbLeft = sidx;
        } else {
            d.numSucceedingNeighbor += numNeighbors;
            d.lastProbRight = sidx;
        }
        return numNeighbors;

    }

    public void checkCoreCountPoints() {
        int countOutlier = 0;
        for (CorePoint c : all_distinct_cores) {
            if (c.isCoveredAllSlides() //                    && (c.getTotalHalfRPoints()+c.getTotal32RPoints()+c.getTotalRPoints() < Constants.k)
                    ) {
                countOutlier += 1;
            }
        }
        System.out.println("Num Core Contains all slides = " + countOutlier);
    }

    public boolean check_distance_neighbor_boolean(C_Data d, C_Data d2) {

        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d <= Constants.R;
    }

    public double check_distance_neighbor(C_Data d, C_Data d2) {

        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    private void processExpiredData(int expiredSlideIndex) {
        all_slides.remove(expiredSlideIndex);

        all_core_points.remove(expiredSlideIndex);
        for (CorePoint c : all_distinct_cores) {
            if (c.closeNeighbors_halfR.get(expiredSlideIndex) != null) {
                c.totalHalfRPoints -= c.closeNeighbors_halfR.get(expiredSlideIndex).size();
                c.closeNeighbors_halfR.remove(expiredSlideIndex);
            }
            c.closeNeighbors_R.remove(expiredSlideIndex);
            c.closeNeighbors_3halfR.remove(expiredSlideIndex);
            c.closeNeighbors_2R.remove(expiredSlideIndex);
        }

    }

    private HashMap<Integer, ArrayList<C_Data>> indexByAtt(ArrayList<C_Data> datas, int i, double min_value) {

        HashMap<Integer, ArrayList<C_Data>> results = new HashMap<>();

        for (C_Data d : datas) {
            int bin = (int) ((d.values[i] - min_value) / (Constants.R / 2));
            if (results.containsKey(bin)) {
                results.get(bin).add(d);
            } else {
                ArrayList<C_Data> v = new ArrayList<>();
                v.add(d);
                results.put(bin, v);
            }
        }
        return results;
    }

    private void scanForCore(CorePoint c, int sIdx,
            HashMap<Integer, ArrayList<C_Data>> bin_map, int index_att, double min_value) {

        ArrayList<C_Data> neighborsInHalfR = new ArrayList<>();
        ArrayList<C_Data> neighborsInR = new ArrayList<>();
        ArrayList<C_Data> neighborsIn3HalfR = new ArrayList<>();
        ArrayList<C_Data> neighborsIn2R = new ArrayList<>();

        int bin_core = (int) ((c.values[index_att] - min_value) / (Constants.R / 2));
        if (bin_map.containsKey(bin_core)) {
            neighborsInHalfR.addAll(bin_map.get(bin_core));
        }
        int[] possible_bins = new int[]{bin_core - 4, bin_core - 3, bin_core - 2, bin_core - 1,
            bin_core + 1, bin_core + 2, bin_core + 3, bin_core + 4};

        for (int b : possible_bins) {
            if (bin_map.containsKey(b)) {
                for (C_Data d2 : bin_map.get(b)) {
                    double distance = DistanceFunction.euclideanDistance(c, d2);
                    if (distance <= Constants.R / 2) {
                        neighborsInHalfR.add(d2);
                        d2.closeCoreMaps_halfR = c;
                    } else if (distance <= Constants.R) {
                        neighborsInR.add(d2);
                        d2.closeCoreMaps_R.add(c);
                    } else if (distance <= Constants.R * 1.5) {
                        neighborsIn3HalfR.add(d2);

                    } else if (distance <= Constants.R * 2) {
                        neighborsIn2R.add(d2);

                    }

                }
            }
        }
//        }
        c.closeNeighbors_halfR.put(sIdx, neighborsInHalfR);
        c.closeNeighbors_R.put(sIdx, neighborsInR);
        c.closeNeighbors_3halfR.put(sIdx, neighborsIn3HalfR);
        c.closeNeighbors_2R.put(sIdx, neighborsIn2R);
    }

    private void scanForCore(CorePoint c, int sIdx) {

        ArrayList<C_Data> neighborsInHalfR = new ArrayList<>();
        ArrayList<C_Data> neighborsInR = new ArrayList<>();
        ArrayList<C_Data> neighborsIn3HalfR = new ArrayList<>();
        ArrayList<C_Data> neighborsIn2R = new ArrayList<>();

        for (C_Data d2 : all_slides.get(sIdx)) {
            double distance = DistanceFunction.euclideanDistance(c, d2);
            if (distance <= Constants.R / 2) {
                neighborsInHalfR.add(d2);
                d2.closeCoreMaps_halfR = c;
            } else if (distance <= Constants.R) {
                neighborsInR.add(d2);
                d2.closeCoreMaps_R.add(c);
            } else if (distance <= Constants.R * 1.5) {
                neighborsIn3HalfR.add(d2);

            } else if (distance <= Constants.R * 2) {
                neighborsIn2R.add(d2);

            }

        }
//        }

        c.closeNeighbors_halfR.put(sIdx, neighborsInHalfR);
        c.closeNeighbors_R.put(sIdx, neighborsInR);
        c.closeNeighbors_3halfR.put(sIdx, neighborsIn3HalfR);
        c.closeNeighbors_2R.put(sIdx, neighborsIn2R);
    }

    private ArrayList<CorePoint> selectCore(Integer sIdx) {

        ArrayList<CorePoint> corePoints = new ArrayList<>();

        for (int i = 0; i < Constants.slide; i++) {

            C_Data d = all_slides.get(sIdx).get(i);

            //using mtree
            MTreeCorePoint.Query query = mtree.getNearestByRange(d, Constants.R * 2);
            for (MTreeClass.ResultItem ri : query) {
                CorePoint c = (CorePoint) ri.data;
                double distance = ri.distance;

                if (distance <= Constants.R / 2) {
                    d.closeCoreMaps_halfR = c;
                    ArrayList<C_Data> closeNeighbors = c.closeNeighbors_halfR.get(sIdx);
                    if (closeNeighbors == null) {
                        closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                        closeNeighbors.add(d);
                        c.closeNeighbors_halfR.put(sIdx, closeNeighbors);
                    } else {
                        closeNeighbors.add(d);
                    }

                } else if (distance <= Constants.R) {

                    d.closeCoreMaps_R.add(c);

                    ArrayList<C_Data> closeNeighbors = c.closeNeighbors_R.get(sIdx);
                    if (closeNeighbors == null) {
                        closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                        closeNeighbors.add(d);
                        c.closeNeighbors_R.put(sIdx, closeNeighbors);
                    } else {
                        closeNeighbors.add(d);
                    }
                } else if (distance <= Constants.R * 2) {
                    d.closeCoreMaps_2R.add(c);

                }
                if (!corePoints.contains(c)) {
                    corePoints.add(c);
                }
//                mtree.add(c);
            }

            //create a new core from d
            if (d.closeCoreMaps_R.isEmpty() && d.closeCoreMaps_halfR == null) {

                CorePoint c = new CorePoint(d);
                all_distinct_cores.add(c);
                mtree.add(c);
                corePoints.add(c);

                ArrayList<C_Data> closeNeighbors = c.closeNeighbors_halfR.get(sIdx);
                if (closeNeighbors == null) {
                    closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                        closeNeighbors.add(d);
                    c.closeNeighbors_halfR.put(sIdx, closeNeighbors);
                } else {
                    closeNeighbors.add(d);
                }
                d.closeCoreMaps_halfR = c;
            }

        }
//        System.out.println("Finished round 1!");
//        System.out.println("Num core points = " + corePoints.size());
        // scan points for cores
        for (CorePoint c : corePoints) {
            if (c.closeNeighbors_halfR.get(sIdx) == null) {
                c.closeNeighbors_halfR.put(sIdx, new ArrayList<>());
            }
            if (c.closeNeighbors_R.get(sIdx) == null) {
                c.closeNeighbors_R.put(sIdx, new ArrayList<>());
            }
            if (c.closeNeighbors_3halfR.get(sIdx) == null) {
                c.closeNeighbors_3halfR.put(sIdx, new ArrayList<>());
            }
            if (c.closeNeighbors_2R.get(sIdx) == null) {
                c.closeNeighbors_2R.put(sIdx, new ArrayList<>());
            }
            for (CorePoint c2 : corePoints) {

                if (c != c2) {
                    double distance = DistanceFunction.euclideanDistance(c, c2);
//                    if (MTTest.numberWindows > 1) {
//                        numDCSForIndexing += 1;
//                    }
                    if (distance <= Constants.R * 3) {
                        probCoreWithList(c, c2.closeNeighbors_halfR.get(sIdx), sIdx);
                        probCoreWithList(c, c2.closeNeighbors_R.get(sIdx), sIdx);
                    }
                }
            }

        }
//        for (CorePoint c : newCores) {
//            mtree.add(c);
//        }
        return corePoints;

    }

    private int probe_current_slide(C_Data d) {
        int count = 0;
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();
        if (d.closeCoreMaps_halfR != null) {
            //prob with candidate 
            count = d.closeCoreMaps_halfR.closeNeighbors_R.get(d.sIndex).size();
//            possibleCandidates.add(d.closeCoreMaps_halfR.closeNeighbors_halfR.get(d.sIndex));
            possibleCandidates.add(d.closeCoreMaps_halfR.closeNeighbors_R.get(d.sIndex));
            possibleCandidates.add(d.closeCoreMaps_halfR.closeNeighbors_3halfR.get(d.sIndex));

        } else {
            CorePoint core = d.closeCoreMaps_R.get(0);
            possibleCandidates.add(core.closeNeighbors_halfR.get(d.sIndex));
            possibleCandidates.add(core.closeNeighbors_R.get(d.sIndex));
            possibleCandidates.add(core.closeNeighbors_3halfR.get(d.sIndex));
            possibleCandidates.add(core.closeNeighbors_2R.get(d.sIndex));
        }

        for (ArrayList<C_Data> ps : possibleCandidates) {
            for (int t = 0; t < ps.size(); t++) {
                C_Data d2 = ps.get(t);

                if (check_distance_neighbor_boolean(d, d2)) {

                    count += 1;
                }
//                   

            }

        }
        d.lastProbRight = d.sIndex;
        d.lastProbLeft = d.sIndex;
        return count;
    }

    private void probCoreWithList(CorePoint c, ArrayList<C_Data> candidates, int sIdx) {
        if (candidates != null) {

            for (C_Data d2 : candidates) {
                double distance = DistanceFunction.euclideanDistance(c, d2);
//                if(MTTest.numberWindows > 1){
//                    numDCSForIndexing +=1;
//                }
                if (distance <= Constants.R / 2) {
                    c.closeNeighbors_halfR.get(sIdx).add(d2);
                    d2.closeCoreMaps_halfR = c;
                } else if (distance <= Constants.R) {
                    c.closeNeighbors_R.get(sIdx).add(d2);
                    d2.closeCoreMaps_R.add(c);
                } else if (distance <= Constants.R * 1.5) {
                    c.closeNeighbors_3halfR.get(sIdx).add(d2);

                } else if (distance <= Constants.R * 2) {
                    c.closeNeighbors_2R.get(sIdx).add(d2);

                }

            }
//        }

//            c.closeNeighbors_halfR.put(sIdx, neighborsInHalfR);
//            c.closeNeighbors_R.put(sIdx, neighborsInR);
//            c.closeNeighbors_3halfR.put(sIdx, neighborsIn3HalfR);
//            c.closeNeighbors_2R.put(sIdx, neighborsIn2R);
        }
    }

    final class ResultFindCore {

        private final double distance;
        private final ArrayList<CorePoint> cores;
        public ArrayList<Double> distance_to_cores = new ArrayList<>();

        public ResultFindCore(double distance, ArrayList<CorePoint> cores) {
            this.distance = distance;
            this.cores = cores;
        }

        public ResultFindCore(double distance, ArrayList<CorePoint> cores, ArrayList<Double> all_distances) {
            this.distance = distance;
            this.cores = cores;
            this.distance_to_cores = all_distances;
        }

        public double getDistance() {
            return this.distance;
        }

        public ArrayList<CorePoint> getCore() {
            return this.cores;
        }
    }

    private ResultFindCore findCloseCore(C_Data d, int slideIndex) {

        ArrayList<CorePoint> resultCore = null;

        if (d.closeCoreMaps_halfR != null
                && d.closeCoreMaps_halfR.closeNeighbors_halfR.containsKey(slideIndex)) {
            resultCore = new ArrayList<>();
            resultCore.add(d.closeCoreMaps_halfR);
            return new ResultFindCore(Constants.R / 2, resultCore);
        } else if (!d.closeCoreMaps_R.isEmpty()) {
            for (CorePoint c : d.closeCoreMaps_R) {
                if (c.closeNeighbors_2R.containsKey(slideIndex)) {
                    resultCore = new ArrayList<>();
                    resultCore.add(c);
                    return new ResultFindCore(Constants.R, resultCore);
                }
            }

        }

        //range query to find close core in range R/2
//        MTreeCorePoint.Query query = all_indexed_cores.get(slideIndex).getNearestByRange(d, Constants.R / 2);
//        ArrayList<CorePoint> queryResult = new ArrayList<>();
//        for (MTreeCorePoint.ResultItem ri : query) {
//            queryResult.add((CorePoint) ri.data);
//        }
//        if (!queryResult.isEmpty()) {
//            return new ResultFindCore(Constants.R / 2, queryResult);
//        } else {
//            //range query in range R
//            query = all_indexed_cores.get(slideIndex).getNearestByRange(d, Constants.R);
//
//            for (MTreeCorePoint.ResultItem ri : query) {
//                queryResult.add((CorePoint) ri.data);
//            }
//            if (!queryResult.isEmpty()) {
//                return new ResultFindCore(Constants.R, queryResult);
//            } else {
//                query = all_indexed_cores.get(slideIndex).getNearestByRange(d, Constants.R * 2);
//                for (MTreeCorePoint.ResultItem ri : query) {
//                    queryResult.add((CorePoint) ri.data);
//                }
//                if (!queryResult.isEmpty()) {
//                    return new ResultFindCore(Constants.R * 2, queryResult);
//                } else {
//                    return new ResultFindCore(Constants.R * 2, null);
//                }
//            }
//        }
        ArrayList<CorePoint> corePoints = all_core_points.get(slideIndex);

        ArrayList<CorePoint> inRangeRCores = new ArrayList<>();
        ArrayList<CorePoint> inRangeDoubleRCores = new ArrayList<>();
        ArrayList<Double> distance_to_cores = new ArrayList<>();
        if (corePoints != null) {
            for (int i = 0; i < corePoints.size(); i++) {
                CorePoint c = corePoints.get(i);
                double distance = DistanceFunction.euclideanDistance(d, c);
//                if (count > 1) {
//                    numDCS += 1;
//                }
//            if (distance <= Constants.R / 2) {
//                resultCore = new ArrayList<>();
//                resultCore.add(c);
//                return new ResultFindCore(Constants.R / 2, resultCore);
//            } else 
                if (distance <= Constants.R) {
                    inRangeRCores.add(c);
                    break;
                    //test
//                return new ResultFindCore(Constants.R, inRangeRCores);
                    //end test
                } else if (distance <= Constants.R * 2) {
                    inRangeDoubleRCores.add(c);
                    distance_to_cores.add(distance);
                }
            }
        }
        if (!inRangeRCores.isEmpty()) {
//            System.out.println("in Range R core = "+ inRangeRCores.size());

            return new ResultFindCore(Constants.R, inRangeRCores);
        } else if (!inRangeDoubleRCores.isEmpty()) {
//            System.out.println("AAAAAAAAAAAAAa");
            return new ResultFindCore(Constants.R * 2, inRangeDoubleRCores, distance_to_cores);
        } else {
            return new ResultFindCore(Constants.R * 2, null);

        }
    }

    class C_Data extends Data {

        private int numSucceedingNeighbor = 0;
        private int numNeighborSameSlide = 0;
        public int lastProbRight = -1;
        public int lastProbLeft = -1;
//        private double supportProb = 0;
        public int neighborCount = 0;
        public boolean isSafe = false;
        public CorePoint closeCoreMaps_halfR;
        public ArrayList<CorePoint> closeCoreMaps_R = new ArrayList<>();
        public ArrayList<CorePoint> closeCoreMaps_2R = new ArrayList<>();
        public int sIndex = -1;

        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;

            this.sIndex = (arrivalTime - 1) / Constants.slide;

        }

        public C_Data() {

        }

        public int countNeighbor() {
            return neighborCount;
        }

    }

    class CorePoint extends C_Data implements Comparable<Data> {

        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_halfR = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_R = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_3halfR = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_2R = new HashMap<>();

        public int totalHalfRPoints = 0;
//        public int totalRPoints = 0;
        public int total2RPoints = 0;
//        public int total3HalfRPoints = 0;

        public int getTotalHalfRPoints() {
            int t = 0;
            t = closeNeighbors_halfR.entrySet().parallelStream().map((e) -> e.getValue().size()).reduce(t, Integer::sum);
            return t;
        }

        public int getTotal32RPoints() {
            int t = 0;
            for (Map.Entry<Integer, ArrayList<C_Data>> e : closeNeighbors_3halfR.entrySet()) {
                t += e.getValue().size();
            }
            return t;
        }

        public int getTotal2RPoints() {
            int t = 0;
            for (Map.Entry<Integer, ArrayList<C_Data>> e : closeNeighbors_2R.entrySet()) {
                t += e.getValue().size();
            }
            return t;
        }

        public int getTotalRPoints() {
            int t = 0;
            for (Map.Entry<Integer, ArrayList<C_Data>> e : closeNeighbors_2R.entrySet()) {
                t += e.getValue().size();
            }
            return t;
        }

        public boolean isCoveredAllSlides() {
            for (Map.Entry<Integer, ArrayList<C_Data>> e : all_slides.entrySet()) {
                if (!closeNeighbors_halfR.containsKey(e.getKey())) {
                    return false;
                }
            }
            return true;
        }

//        public int getTotalCoverPoint() {
//            return closeNeighbors_2R.size() + closeNeighbors_3halfR.size()
//                    + closeNeighbors_R.size() + closeNeighbors_halfR.size();
//        }
//  
//        private int numPointsCovered = 0;
        public CorePoint(C_Data d) {
            this.values = d.values;
            this.hashCode = d.hashCode;
            this.arrivalTime = d.arrivalTime;
        }

    }

}
