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
import java.util.Set;
import mtree.ComposedSplitFunction;
import mtree.DistanceFunctions;
import mtree.MTree;
import mtree.PartitionFunctions;
import mtree.PromotionFunction;
import mtree.tests.Data;
import mtree.tests.MTTest;
import mtree.utils.Constants;
import mtree.utils.Pair;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class Approx_CPOD {

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

//    public static double timeProcessExpiredSlide = 0;
//    public static double timeCreatingCore = 0;
//    public static double timeProbing = 0;
//    public static double timeReProbing = 0;
//    public static int count = 0;
//    public static double numPointNeedProb = 0;
//    public static double avg_points_check = 0;
//    public static int countPoint = 0;
//    public static long numPointNeedNS = 0;
//    public static long numDCS = 0;
//    public static long count = 0;
//    public static long numDCSForIndexing = 0;
//    public static int maxMembers = 15;
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
//            count += 1;
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

        long start = Utils.getCPUTime();
        ArrayList<Data> result = new ArrayList<>((int) (data.size() * 1.0 / 100));
        expiredSlideIndex = (currentTime - 1) / slide - Constants.W / slide;
//        System.out.println("Expire slide index = " + expiredSlideIndex);
//        long start = Utils.getCPUTime();
        processExpiredData(expiredSlideIndex);
//        if (data.size() != Constants.W) {
//            count += 1;
//            timeProcessExpiredSlide += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
//        }
//        System.out.println("Time for processing expired slide = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        start = Utils.getCPUTime();
        for (int sIdx : slide_to_process) {

//            HashMap<Integer, ArrayList<C_Data>> bin_map = indexByAtt(all_slides.get(sIdx), 0, 0);
            ArrayList<CorePoint> corePoints = selectCore(sIdx);
//            for (CorePoint c : corePoints) {
//                scanForCore(c, sIdx);
//            }

            // add core points to mtree
//            MTreeCorePoint mtree = new MTreeCorePoint();
//
//            for (CorePoint c : corePoints) {
//                mtree.add(c);
//            }
            all_core_points.put(sIdx, corePoints);
//            all_indexed_cores.put(sIdx, mtree);
//            System.out.println("Num new core = " + corePoints.size());
//            System.out.println("All distinc core = " + all_distinct_cores.size());
        }

//        System.out.println("Time for creating cores = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
//        
//        System.out.println("Probing");
//        start = Utils.getCPUTime();
        int newestSlide = (currentTime - 1) / Constants.slide;
//        timeCheckingCandidates = 0;
//        timeFindingCore = 0;
//        timeAddingToNeighborCount = 0;
//        timeMarkProbing = 0;
//        timeForAddingToOutlierList = 0;

//        checkCoreCountPoints();
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
//        if (count > 0) {
//
//            timeCreatingCore += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
//        }
//        for (CorePoint c : all_distinct_cores) {
//            if (data.size() == Constants.W) {
//                c.totalHalfRPoints = c.getTotalHalfRPoints();
//            } else {
//                if (c.closeNeighbors_halfR.get(newestSlide) != null) {
//                    c.totalHalfRPoints += c.closeNeighbors_halfR.get(newestSlide).size();
//                }
//            }
//
//        }
        if (MTTest.numberWindows > 1) {
            MTTest.timeForIndexing += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        }
//
        start = Utils.getCPUTime();
        int pointNeedProb = 0;

        for (int i = 0; i < d_to_process.size(); i++) {
            C_Data d = d_to_process.get(i);
            if (d.closeCoreMaps_halfR != null && d.closeCoreMaps_halfR.totalHalfRPoints >= Constants.k + 1) {
                continue;
            }
            if (d.neighborCount < Constants.k) {

                probe(d, newestSlide);
                pointNeedProb += 1;
            }
//            System.out.println("Finished "+ i);
        }
//        if (count >= 1) {
//            numPointNeedProb += pointNeedProb;
//        }

//        System.out.println("Points need Prob = " + pointNeedProb);
//        System.out.println("Time finding cores = " + timeFindingCore);
//        System.out.println("Time checking candidates = " + timeCheckingCandidates);
//        System.out.println("Time adding to neighbor counts = " + timeAddingToNeighborCount);
//        System.out.println("Time mark probing = " + timeMarkProbing);
//        System.out.println("Time for adding to outlier list = " + timeForAddingToOutlierList);
//        System.out.println("Time for probing = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
//        System.out.println("re-probing + outlier report");
//        if (count >= 1) {
//            timeProbing += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
//        }
        start = Utils.getCPUTime();
        for (Map.Entry<Integer, HashSet<C_Data>> e : outlierList.entrySet()) {
//        for (int slideIndex : outlierList.keySet()) {
            for (C_Data d : e.getValue()) {

                if (d.closeCoreMaps_halfR != null && d.closeCoreMaps_halfR.totalHalfRPoints >= Constants.k + 1) {
                    continue;
                }
                if (d.neighborCount < Constants.k && d.sIndex < newestSlide) {
                    if (d.lastProbRight < newestSlide) {

                        probe(d, newestSlide);
                    }

                }
                if (d.neighborCount < Constants.k) {
                    result.add(d);
                }
            }
        }

        if (MTTest.numberWindows > 1) {
            MTTest.timeForNeighborSearch += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        }
//        if (count >= 1) {
//            timeReProbing += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
//        };
//        System.out.println("Time for reprobing+report outliers = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        return result;

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
        outlierList.remove(expiredSlideIndex);
        if (neighborCountTrigger.containsKey(expiredSlideIndex)) {
            for (C_Data d : neighborCountTrigger.get(expiredSlideIndex)) {
                d.neighborCount -= d.pred_neighbor_count.get(expiredSlideIndex);
                d.pred_neighbor_count.remove(expiredSlideIndex);
            }
        }
        if (neighborCountTrigger.containsKey(expiredSlideIndex - 1)) {
            neighborCountTrigger.remove(expiredSlideIndex - 1);
        }

        all_core_points.remove(expiredSlideIndex);
//        all_indexed_cores.remove(expiredSlideIndex);
        for (CorePoint c : all_distinct_cores) {
            if (c.closeNeighbors_halfR.get(expiredSlideIndex) != null) {
                c.totalHalfRPoints -= c.closeNeighbors_halfR.get(expiredSlideIndex).size();
                c.closeNeighbors_halfR.remove(expiredSlideIndex);
            }
            c.closeNeighbors_R.remove(expiredSlideIndex);
            c.countNeighbors_3halfR.remove(expiredSlideIndex);
            c.countNeighbors_halfR.remove(expiredSlideIndex);
            c.countNeighbors_R.remove(expiredSlideIndex);
            c.countNeighbors_2R.remove(expiredSlideIndex);
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

//    private void scanForCore(CorePoint c, int sIdx,
//            HashMap<Integer, ArrayList<C_Data>> bin_map, int index_att, double min_value) {
//
//        ArrayList<C_Data> neighborsInHalfR = new ArrayList<>();
//        ArrayList<C_Data> neighborsInR = new ArrayList<>();
//        ArrayList<C_Data> neighborsIn3HalfR = new ArrayList<>();
//        ArrayList<C_Data> neighborsIn2R = new ArrayList<>();
//
//        int bin_core = (int) ((c.values[index_att] - min_value) / (Constants.R / 2));
//        if (bin_map.containsKey(bin_core)) {
//            neighborsInHalfR.addAll(bin_map.get(bin_core));
//        }
//        int[] possible_bins = new int[]{bin_core - 4, bin_core - 3, bin_core - 2, bin_core - 1,
//            bin_core + 1, bin_core + 2, bin_core + 3, bin_core + 4};
//
//        for (int b : possible_bins) {
//            if (bin_map.containsKey(b)) {
//                for (C_Data d2 : bin_map.get(b)) {
//                    double distance = DistanceFunction.euclideanDistance(c, d2);
//                    if (distance <= Constants.R / 2) {
//                        neighborsInHalfR.add(d2);
//                        d2.closeCoreMaps_halfR = c;
//                    } else if (distance <= Constants.R) {
//                        neighborsInR.add(d2);
//                        d2.closeCoreMaps_R.add(c);
//                    } else if (distance <= Constants.R * 1.5) {
//                        neighborsIn3HalfR.add(d2);
//
//                    } else if (distance <= Constants.R * 2) {
//                        neighborsIn2R.add(d2);
//
//                    }
//
//                }
//            }
//        }
////        }
//        c.closeNeighbors_halfR.put(sIdx, neighborsInHalfR);
//        c.closeNeighbors_R.put(sIdx, neighborsInR);
//        c.closeNeighbors_3halfR.put(sIdx, neighborsIn3HalfR);
//        c.closeNeighbors_2R.put(sIdx, neighborsIn2R);
//    }
//    private void scanForCore(CorePoint c, int sIdx) {
//
//        ArrayList<C_Data> neighborsInHalfR = new ArrayList<>();
//        ArrayList<C_Data> neighborsInR = new ArrayList<>();
//        ArrayList<C_Data> neighborsIn3HalfR = new ArrayList<>();
//        ArrayList<C_Data> neighborsIn2R = new ArrayList<>();
//
//        for (C_Data d2 : all_slides.get(sIdx)) {
//            double distance = DistanceFunction.euclideanDistance(c, d2);
//            if (distance <= Constants.R / 2) {
//                neighborsInHalfR.add(d2);
//                d2.closeCoreMaps_halfR = c;
//            } else if (distance <= Constants.R) {
//                neighborsInR.add(d2);
//                d2.closeCoreMaps_R.add(c);
//            } else if (distance <= Constants.R * 1.5) {
//                neighborsIn3HalfR.add(d2);
//
//            } else if (distance <= Constants.R * 2) {
//                neighborsIn2R.add(d2);
//
//            }
//
//        }
////        }
//
//        c.closeNeighbors_halfR.put(sIdx, neighborsInHalfR);
//        c.closeNeighbors_R.put(sIdx, neighborsInR);
//        c.closeNeighbors_3halfR.put(sIdx, neighborsIn3HalfR);
//        c.closeNeighbors_2R.put(sIdx, neighborsIn2R);
//    }
    private ArrayList<CorePoint> selectCore(Integer sIdx) {

        ArrayList<CorePoint> corePoints = new ArrayList<>();
        ArrayList<CorePoint> newCores = new ArrayList<>();

//        ArrayList<C_Data> temp_list = new ArrayList<>();
//        MTreeCorePoint mt = new MTreeCorePoint();
        for (int i = 0; i < Constants.slide; i++) {
//            System.out.println("Slide index = "+ sIdx);

            C_Data d = all_slides.get(sIdx).get(i);

            //scan with current cores first
            for (int j = corePoints.size() - 1; j >= 0; j--) {
                CorePoint c = corePoints.get(j);
                double distance = DistanceFunction.euclideanDistance(d, c);
//                if (MTTest.numberWindows > 1) {
//                    numDCSForIndexing += 1;
//                }

                if (distance <= Constants.R / 2) {
                    ArrayList<C_Data> closeNeighbors = c.closeNeighbors_halfR.get(sIdx);
                    if (closeNeighbors == null) {
                        closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                        closeNeighbors.add(d);
                        c.closeNeighbors_halfR.put(sIdx, closeNeighbors);
                    } else {
                        closeNeighbors.add(d);
                    }
                    d.closeCoreMaps_halfR = c;
                    break;
                } else if (distance <= Constants.R) {
                    ArrayList<C_Data> closeNeighbors = c.closeNeighbors_R.get(sIdx);
                    if (closeNeighbors == null) {
                        closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                        closeNeighbors.add(d);
                        c.closeNeighbors_R.put(sIdx, closeNeighbors);
                    } else {
                        closeNeighbors.add(d);
                    }
                    d.closeCoreMaps_R.add(c);
                    break;
                }
            }

            if ((d.closeCoreMaps_R.isEmpty() && d.closeCoreMaps_halfR == null)) {

//                ArrayList<CorePoint> inRCore = new ArrayList<>();
                //using mtree
                MTreeCorePoint.Query query = mtree.getNearest(d, Constants.R, 1);
                CorePoint c = null;
                double distance = Double.MAX_VALUE;
                for (MTreeClass.ResultItem ri : query) {
                    c = (CorePoint) ri.data;
                    distance = ri.distance;
                }
                if (distance <= Constants.R) {
                    //add c to the core of this slide
                    corePoints.add(c);

                    if (distance <= Constants.R / 2) {
                        ArrayList<C_Data> closeNeighbors = c.closeNeighbors_halfR.get(sIdx);
                        if (closeNeighbors == null) {
                            closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                            closeNeighbors.add(d);
                            c.closeNeighbors_halfR.put(sIdx, closeNeighbors);
                        } else {
                            closeNeighbors.add(d);
                        }
                        d.closeCoreMaps_halfR = c;

                    } else {
                        ArrayList<C_Data> closeNeighbors = c.closeNeighbors_R.get(sIdx);
                        if (closeNeighbors == null) {
                            closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                            closeNeighbors.add(d);
                            c.closeNeighbors_R.put(sIdx, closeNeighbors);
                        } else {
                            closeNeighbors.add(d);
                        }
                        d.closeCoreMaps_R.add(c);

                    }
//                    scanForCore(c, sIdx);

                } else {
                    c = new CorePoint(d);
                    all_distinct_cores.add(c);
                    newCores.add(c);
//                    mtree.add(c);
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

                    //probe neighbors for c
//                    scanForCore(c, sIdx);
                }

            }
        }

        //find scan for cores
        for (CorePoint c : corePoints) {
            if (c.closeNeighbors_halfR.get(sIdx) == null) {
                c.closeNeighbors_halfR.put(sIdx, new ArrayList<>());
            }
            if (c.closeNeighbors_R.get(sIdx) == null) {
                c.closeNeighbors_R.put(sIdx, new ArrayList<>());
            }
            if (c.countNeighbors_R.get(sIdx) == null) {
                c.countNeighbors_R.put(sIdx, new CompactRange(0, 0));
            }
            if (c.countNeighbors_halfR.get(sIdx) == null) {
                c.countNeighbors_halfR.put(sIdx, new CompactRange(0, 0));
            }
            if (c.countNeighbors_3halfR.get(sIdx) == null) {
                c.countNeighbors_3halfR.put(sIdx, new CompactRange(0, 0));
            }
            if (c.countNeighbors_2R.get(sIdx) == null) {
                c.countNeighbors_2R.put(sIdx, new CompactRange(0, 0));
            }

//            c.totalPoints = c.getTotalPoints();
//
//            if (c.totalPoints < maxMembers) {
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
//            }

        }
        for (CorePoint c : newCores) {
            mtree.add(c);
        }
        return corePoints;

    }

    private void probe_slide_right(C_Data d, int slideIndex) {
//        int countNeighbor = 0;

        //find close core
//        long start = Utils.getCPUTime();
        ResultFindCore rf = findCloseCore(d, slideIndex);
//        System.out.println("Time to find core = "+ (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        int case_ = 0;
        if (cores != null) {
            if (distance <= Constants.R / 2) {
                CorePoint c = cores.get(0);

                int possibleNeighbors = c.countNeighbors_halfR.get(slideIndex).count + c.countNeighbors_R.get(slideIndex).count
                        + c.countNeighbors_3halfR.get(slideIndex).count;

                int estimatedCountNeighbor = (int) (possibleNeighbors
                        / Math.pow(c.countNeighbors_3halfR.get(slideIndex).max_dist / Constants.R, c.dimensions()));
//                System.out.println("Estimated Count neighbor = "+ estimatedCountNeighbor);
                if (estimatedCountNeighbor > c.countNeighbors_halfR.get(slideIndex).count) {
                    d.neighborCount += estimatedCountNeighbor;
                    d.numSucceedingNeighbor += estimatedCountNeighbor;
                } else {
                    d.neighborCount += c.countNeighbors_halfR.get(slideIndex).count;
                    d.numSucceedingNeighbor += c.countNeighbors_halfR.get(slideIndex).count;
                }

            } else if (distance <= Constants.R) {
                CorePoint c = cores.get(0);
                int possibleNeighbors = c.countNeighbors_halfR.get(slideIndex).count + c.countNeighbors_R.get(slideIndex).count
                        + c.countNeighbors_3halfR.get(slideIndex).count + c.countNeighbors_2R.get(slideIndex).count;
                int estimatedCountNeighbor = (int) (possibleNeighbors
                        / Math.pow(c.countNeighbors_2R.get(slideIndex).max_dist / Constants.R, c.dimensions()));
                d.neighborCount += estimatedCountNeighbor;
                d.numSucceedingNeighbor += estimatedCountNeighbor;

            } else if (distance <= Constants.R * 2) {

                int estimatedCountNeighbor = 0;
                for (CorePoint c : cores) {
                    int possibleNeighbors = c.countNeighbors_halfR.get(slideIndex).count + c.countNeighbors_R.get(slideIndex).count
                            + c.countNeighbors_3halfR.get(slideIndex).count + c.countNeighbors_2R.get(slideIndex).count;
                    estimatedCountNeighbor += (int) (possibleNeighbors / 2
                            / Math.pow(c.countNeighbors_2R.get(slideIndex).max_dist / Constants.R, c.dimensions()));

                    d.neighborCount += estimatedCountNeighbor;
                    d.numSucceedingNeighbor += estimatedCountNeighbor;
                    //grab close neighbor in range R/2 of c
                    if (d.numSucceedingNeighbor >= Constants.k) {

                        return;

                    }
                }
            }

        }
    }

    private void probe_slide_left(C_Data d, int slideIndex) {

//        int countNeighbor = 0;
        int oldNumNeighbor = d.neighborCount;
        //scan possible points 
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();

        //find close core
//        long start = Utils.getCPUTime();
        ResultFindCore rf = findCloseCore(d, slideIndex);
//        System.out.println("Time to find core = "+ (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        int case_ = 0;
        if (cores != null) {
            if (distance <= Constants.R / 2) {
                CorePoint c = cores.get(0);

                int possibleNeighbors = c.countNeighbors_halfR.get(slideIndex).count + c.countNeighbors_R.get(slideIndex).count
                        + c.countNeighbors_3halfR.get(slideIndex).count;

                int estimatedCountNeighbor = (int) (possibleNeighbors
                        / Math.pow(c.countNeighbors_3halfR.get(slideIndex).max_dist / Constants.R, c.dimensions()));

                if (estimatedCountNeighbor > c.countNeighbors_halfR.get(slideIndex).count) {
                    d.neighborCount += estimatedCountNeighbor;
                    d.pred_neighbor_count.put(slideIndex, estimatedCountNeighbor);
                } else {
                    d.neighborCount += c.countNeighbors_halfR.get(slideIndex).count;
                    d.pred_neighbor_count.put(slideIndex, c.countNeighbors_halfR.get(slideIndex).count);
                }

                if (neighborCountTrigger.containsKey(slideIndex)) {
                    neighborCountTrigger.get(slideIndex).add(d);
                } else {
                    HashSet<C_Data> hs = new HashSet<>();
                    hs.add(d);
                    neighborCountTrigger.put(slideIndex, hs);
                }

            } else if (distance <= Constants.R) {

                CorePoint c = cores.get(0);
                int possibleNeighbors = c.countNeighbors_halfR.get(slideIndex).count + c.countNeighbors_R.get(slideIndex).count
                        + c.countNeighbors_3halfR.get(slideIndex).count + c.countNeighbors_2R.get(slideIndex).count;
                int estimatedCountNeighbor = (int) (possibleNeighbors
                        / Math.pow(c.countNeighbors_2R.get(slideIndex).max_dist / Constants.R, c.dimensions()));
                d.neighborCount += estimatedCountNeighbor;
                d.pred_neighbor_count.put(slideIndex, estimatedCountNeighbor);
                if (neighborCountTrigger.containsKey(slideIndex)) {
                    neighborCountTrigger.get(slideIndex).add(d);
                } else {
                    HashSet<C_Data> hs = new HashSet<>();
                    hs.add(d);
                    neighborCountTrigger.put(slideIndex, hs);
                }

            } else if (distance <= Constants.R * 2) {
                int estimatedCountNeighbor = 0;
                for (CorePoint c : cores) {
                    int possibleNeighbors = c.countNeighbors_halfR.get(slideIndex).count + c.countNeighbors_R.get(slideIndex).count
                            + c.countNeighbors_3halfR.get(slideIndex).count + c.countNeighbors_2R.get(slideIndex).count;
                    estimatedCountNeighbor += (int) (possibleNeighbors / 2
                            / Math.pow(c.countNeighbors_2R.get(slideIndex).max_dist / Constants.R, c.dimensions()));

                    d.neighborCount += estimatedCountNeighbor;
                    d.pred_neighbor_count.put(slideIndex, estimatedCountNeighbor);
                    if (neighborCountTrigger.containsKey(slideIndex)) {
                        neighborCountTrigger.get(slideIndex).add(d);
                    } else {
                        HashSet<C_Data> hs = new HashSet<>();
                        hs.add(d);
                        neighborCountTrigger.put(slideIndex, hs);
                    }

                    if (d.neighborCount >= Constants.k) {

                        return;

                    }
                }
            }

        }

    }

    private void probe(C_Data d, int newestSlide) {
//        d.neighborCount = d.numSucceedingNeighbor;
        //grab neighbor of close core
//        for (Integer sIdx : all_slides.keySet()) {
//            if (sIdx < d.lastProbLeft) {
//                //collect from core 
//                if (d.closeCoreMaps_halfR != null
//                        && d.closeCoreMaps_halfR.closeNeighbors_halfR.containsKey(sIdx)) {
//                    d.neighborCount += d.closeCoreMaps_halfR.closeNeighbors_halfR.size();
//                    d.pred_neighbor_count.put(sIdx, d.closeCoreMaps_halfR.closeNeighbors_halfR.get(sIdx).size());
//                }
//            } else if (sIdx < d.sIndex && d.pred_neighbor_count.containsKey(sIdx)) {
//                d.neighborCount += d.pred_neighbor_count.get(sIdx);
//
//            } else if (sIdx > d.lastProbRight) {
//                if (d.closeCoreMaps_halfR != null
//                        && d.closeCoreMaps_halfR.closeNeighbors_halfR.containsKey(sIdx)) {
//                    d.neighborCount += d.closeCoreMaps_halfR.closeNeighbors_halfR.get(sIdx).size();
//                }
//            }
//        }
//        if (d.neighborCount < Constants.k) {
//            if (count >= 1) {
//                countPoint += 1;
//            }
//        }
        boolean counted = false;
        if (d.lastProbRight < newestSlide) {
            //prob right first
            int slideIndex = d.lastProbRight + 1;
            if (d.lastProbRight == -1) {
                slideIndex = d.sIndex;
            }
            while (slideIndex <= newestSlide && d.neighborCount < Constants.k) {
//                    if (d.closeCoreMaps_halfR != null
//                            && d.closeCoreMaps_halfR.closeNeighbors_halfR.containsKey(slideIndex)) {
//                        d.neighborCount -= d.closeCoreMaps_halfR.closeNeighbors_halfR.get(slideIndex).size();
//                    }
                if (!counted) {
//                    if (count >= 1) {
//                        numPointNeedNS += 1;
//                    }
                    counted = true;
                }
                probe_slide_right(d, slideIndex);
                d.lastProbRight = slideIndex;
                slideIndex++;
            }
        }
        //prob left
        if (d.neighborCount < Constants.k) {
            int slideIndex = d.lastProbLeft - 1;
            if (d.lastProbLeft == -1) {
                slideIndex = d.sIndex - 1;
            }

            while (slideIndex > expiredSlideIndex && slideIndex >= 0
                    && d.neighborCount < Constants.k) {
//                  while(d.lastProbLeft == -1 || (d.lastProbLeft >= 0 && d.lastProbLeft > expiredSlideIndex))
//                    if (d.closeCoreMaps_halfR != null
//                            && d.closeCoreMaps_halfR.closeNeighbors_halfR.containsKey(slideIndex)) {
//                        d.neighborCount -= d.closeCoreMaps_halfR.closeNeighbors_halfR.get(slideIndex).size();
//                    }
                if (!counted) {
//                    if (count > 1) {
//                        numPointNeedNS += 1;
//                    }
                    counted = true;
                }
                probe_slide_left(d, slideIndex);
                d.lastProbLeft = slideIndex;
                slideIndex--;
            }
        }
//        }
//        long start = Utils.getCPUTime();
        if (d.neighborCount < Constants.k) {
            //add to outlier List
            if (outlierList.containsKey(d.sIndex)) {
                outlierList.get(d.sIndex).add(d);
            } else {
                HashSet hs = new HashSet();
                hs.add(d);
                outlierList.put(d.sIndex, hs);
            }
        }
//        timeForAddingToOutlierList += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
    }

    private void probCoreWithList(CorePoint c, ArrayList<C_Data> candidates, int sIdx) {
        if (candidates != null) {

            for (C_Data d2 : candidates) {

                double distance = DistanceFunction.euclideanDistance(c, d2);
//                if(MTTest.numberWindows > 1){
//                    numDCSForIndexing +=1;
//                }
                if (distance <= Constants.R / 2) {
//                    c.closeNeighbors_halfR.get(sIdx).add(d2);
                    c.countNeighbors_halfR.get(sIdx).count += 1;
                    c.countNeighbors_halfR.get(sIdx).max_dist
                            = Math.max(c.countNeighbors_halfR.get(sIdx).max_dist, distance);
                    d2.closeCoreMaps_halfR = c;

//                    c.totalPoints += 1;
                } else if (distance <= Constants.R) {
                    c.countNeighbors_R.get(sIdx).count += 1;
                    c.countNeighbors_R.get(sIdx).max_dist
                            = Math.max(c.countNeighbors_R.get(sIdx).max_dist, distance);
                    d2.closeCoreMaps_R.add(c);
//                    c.totalPoints += 1;
                } else if (distance <= Constants.R * 1.5) {
                    c.countNeighbors_3halfR.get(sIdx).count += 1;
                    c.countNeighbors_3halfR.get(sIdx).max_dist
                            = Math.max(c.countNeighbors_3halfR.get(sIdx).max_dist, distance);
//                    c.totalPoints += 1;
                } else if (distance <= Constants.R * 2) {
                    c.countNeighbors_2R.get(sIdx).count += 1;
                    c.countNeighbors_2R.get(sIdx).max_dist
                            = Math.max(c.countNeighbors_2R.get(sIdx).max_dist, distance);
//                    c.totalPoints += 1;
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
                if (c.countNeighbors_2R.containsKey(slideIndex)) {
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

//        public boolean addedForProbing = false;
        private int numSucceedingNeighbor = 0;
//        private boolean isOutlier;
        public int lastProbRight = -1;
        public int lastProbLeft = -1;
//        private int lastProbCoreLeft = -1;
//        private int lastProbCoreRight = -1;
        private HashMap<Integer, Integer> pred_neighbor_count = new HashMap<>();
        public int neighborCount = 0;

        private CorePoint closeCoreMaps_halfR;
        private ArrayList<CorePoint> closeCoreMaps_R = new ArrayList<>();
//        private CorePoint closeCore;
//        public ArrayList<CorePoint> linkCoreMaps = new ArrayList<>();

//        public ArrayList<Integer[]> groups = new ArrayList<>();
//        public double[] mean_list;
//        public double[] std_list;
//        public ArrayList<Integer> ordered_attributes = new ArrayList<>();
        public int sIndex = -1;

        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;

            this.sIndex = (arrivalTime - 1) / Constants.slide;

//            if (null != Constants.dataFile) {
//                switch (Constants.dataFile) {
//
//                    case "household2.txt":
//
//                        if (Constants.useLB1 || Constants.useUB1) {
//                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6});
//                        }
//                        if (Constants.useLB2) {
//                            this.ordered_attributes = new ArrayList<>(Arrays.asList(2, 3, 6, 0, 5, 1, 4));
//                        }
//                        break;
//                    case "new_hpc.txt":
//                        this.groups.add(new Integer[]{0, 1, 3});
//                        this.groups.add(new Integer[]{2, 4, 5, 6});
//                        break;
//                    case "covtype.data":
//                        if (Constants.useLB1 || Constants.useUB1) {
////                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
////                                21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
////                                41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54});
////                            this.groups.add(new Integer[]{ });
//                            this.groups.add(new Integer[]{0});
//                            this.groups.add(new Integer[]{1});
//                            this.groups.add(new Integer[]{2});
//                            this.groups.add(new Integer[]{4});
//                            this.groups.add(new Integer[]{3});
//                            this.groups.add(new Integer[]{5});
//                            this.groups.add(new Integer[]{6});
//                            this.groups.add(new Integer[]{7});
//                            this.groups.add(new Integer[]{8});
//                            this.groups.add(new Integer[]{9});
//                            this.groups.add(new Integer[]{11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54});
////                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5});
////                            this.groups.add(new Integer[]{6, 7, 8, 9, 54});
//                        }
//
//                        if (Constants.useLB2) {
//                            this.ordered_attributes = new ArrayList<>(Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54));
//                        }
//                        break;
//                    case "ethylene.txt":
//
//                        if (Constants.useLB1 || Constants.useUB1) {
//
//                            this.groups.add(new Integer[]{0, 2});
////                    this.groups.add(new Integer[]{1});
//                            this.groups.add(new Integer[]{3, 4, 6, 7});
//                            this.groups.add(new Integer[]{8, 14, 15});
//                            this.groups.add(new Integer[]{1, 9});
//                            this.groups.add(new Integer[]{10, 11, 12});
//                            this.groups.add(new Integer[]{13});
//                            this.groups.add(new Integer[]{5});
//                        }
//
//                        if (Constants.useLB2) {
//                            this.ordered_attributes = new ArrayList<>(Arrays.asList(8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7));
//                        }
//                        break;
//
//                    case "new_tao.txt":
//                        this.groups.add(new Integer[]{0, 1, 2});
//
//                        break;
//
//                    case "tao.txt":
//                        if (Constants.useLB1 || Constants.useUB1) {
//                            this.groups.add(new Integer[]{0, 1, 2});
//                        }
//                        if (Constants.useLB2) {
//                            this.ordered_attributes = new ArrayList<>(Arrays.asList(1, 2));
//                        }
//                        break;
//
//                    case "gaussian.txt":
//                        if (Constants.useLB1 || Constants.useUB1) {
//                            this.groups.add(new Integer[0]);
//                        }
//                    default:
//                        break;
//                }
//
//                if (this.groups.size() > 0) {
//                    int num_group = this.groups.size();
//
//                    this.mean_list = new double[num_group];
//                    this.std_list = new double[num_group];
//                    for (int i = 0; i < num_group; i++) {
//                        double temp = 0;
//                        for (Integer j : this.groups.get(i)) {
//                            temp += this.values[j];
//                        }
//                        this.mean_list[i] = temp / this.groups.get(i).length;
//                        if (this.groups.get(i).length == 1) {
//                            this.std_list[i] = 0;
//                        } else {
//                            temp = 0;
//                            for (Integer j : this.groups.get(i)) {
//                                temp += (this.values[j] - this.mean_list[i]) * (this.values[j] - this.mean_list[i]);
//                            }
//
//                            this.std_list[i] = Math.sqrt(temp / this.groups.get(i).length);
//                        }
//                    }
//                }
//            }
//            this.isOutlier = true;
        }

        public C_Data() {

        }

        public int countNeighbor() {
            return neighborCount;
        }

    }

    class CompactRange {

        public int count;
        public double max_dist;

        public CompactRange(int count, double max_dist) {
            this.count = count;
            this.max_dist = max_dist;
        }
    }

    class CorePoint extends C_Data implements Comparable<Data> {

        public HashMap<Integer, CompactRange> countNeighbors_halfR = new HashMap<>();
        public HashMap<Integer, CompactRange> countNeighbors_R = new HashMap<>();
        public HashMap<Integer, CompactRange> countNeighbors_3halfR = new HashMap<>();
        public HashMap<Integer, CompactRange> countNeighbors_2R = new HashMap<>();

        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_halfR = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_R = new HashMap<>();

        public int totalHalfRPoints = 0;
        public int totalRPoints = 0;
        public int total2RPoints = 0;
        public int total3HalfRPoints = 0;
        public int totalPoints = 0;

        public int getTotalPoints() {
            return getTotalHalfRPoints() + getTotal32RPoints() + getTotal2RPoints() + getTotalRPoints();
        }

        public int getTotalHalfRPoints() {
            int t = 0;
            for (Map.Entry<Integer, CompactRange> e : countNeighbors_halfR.entrySet()) {
                t += e.getValue().count;
            }
            return t;
        }

        public int getTotal32RPoints() {
            int t = 0;
            for (Map.Entry<Integer, CompactRange> e : countNeighbors_3halfR.entrySet()) {
                t += e.getValue().count;
            }
            return t;
        }

        public int getTotal2RPoints() {
            int t = 0;
            for (Map.Entry<Integer, CompactRange> e : countNeighbors_2R.entrySet()) {
                t += e.getValue().count;
            }
            return t;
        }

        public int getTotalRPoints() {
            int t = 0;
            for (Map.Entry<Integer, CompactRange> e : countNeighbors_R.entrySet()) {
                t += e.getValue().count;
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
//            this.groups = d.groups;
//            this.mean_list = d.mean_list;
//            this.std_list = d.std_list;
        }

    }

}

//class MTreeCorePoint extends MTree<Data> {
//
//    private static final PromotionFunction<Data> nonRandomPromotion = new PromotionFunction<Data>() {
//        @Override
//        public Pair<Data> process(Set<Data> dataSet, mtree.DistanceFunction<? super Data> distanceFunction) {
//            return mtree.utils.Utils.minMax(dataSet);
//        }
//    };
//
//    MTreeCorePoint() {
//        super(100, DistanceFunctions.EUCLIDEAN, new ComposedSplitFunction<Data>(nonRandomPromotion,
//                new PartitionFunctions.BalancedPartition<Data>()));
//    }
//
//    @Override
//    public void add(Data data) {
//        super.add(data);
//        _check();
//    }
//
//    @Override
//    public boolean remove(Data data) {
//        boolean result = super.remove(data);
//        _check();
//        return result;
//    }
//
//};
