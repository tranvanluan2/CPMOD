/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import mtree.tests.Data;
import mtree.utils.Constants;

/**
 *
 * @author luan
 */
public class CPOD_MQ_ShareCore {

    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public static HashMap<OD_Query, HashMap<Integer, ArrayList<CorePoint>>> all_core_points_map = new HashMap<>();
    public static HashMap<OD_Query, MTreeCorePoint> all_mtree_map = new HashMap<>();

    public static HashMap<OD_Query, ArrayList<CorePoint>> all_distinct_core_map = new HashMap<>();
    public static HashMap<OD_Query, HashMap<Integer, HashSet<C_Data>>> outlierList_map = new HashMap<>();
    public static ArrayList<CorePoint> all_distince_core_point = new ArrayList<>();

    public ArrayList<OD_Query> all_queries = new ArrayList<>();

    public static HashMap<OD_Query, HashMap<Integer, HashSet<C_Data>>> neighborCountTrigger = new HashMap<>();

    public void add_query(OD_Query q) {
        all_queries.add(q);
        
    }

    public void remove_query(OD_Query q) {
        all_queries.remove(q);
    }

    public void remove_query(int idx) {
        all_queries.remove(idx);
    }

    public void update_query(int idx, OD_Query new_q) {
        all_queries.set(idx, new_q);
    }

    public HashMap<OD_Query, HashSet<Data>> slide_process(ArrayList<Data> data, int _currentTime) {

        HashMap<OD_Query, HashSet<Data>> result = new HashMap<>();
        for (OD_Query q : all_queries) {
            result.put(q, new HashSet<>());
        }
        currentTime = _currentTime;
        ArrayList<C_Data> d_to_process = new ArrayList<>(data.size());
        int[] slide_to_process;
        if (currentTime == Constants.W) {
            slide_to_process = new int[Constants.W / Constants.slide];
            for (int i = 0; i < slide_to_process.length; i++) {
                slide_to_process[i] = i;
            }
        } else {
            slide_to_process = new int[]{(currentTime - 1) / Constants.slide};
        }

        //insert new slide
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

        expiredSlideIndex = (currentTime - 1) / Constants.slide - Constants.W / Constants.slide;
        processExpiredData(expiredSlideIndex);

        for (int sIdx : slide_to_process) {
            HashMap<OD_Query, ArrayList<CorePoint>> corepoints = selectCoreForAll(sIdx, all_queries);
            for (OD_Query q : all_queries) {
                if (all_core_points_map.get(q) == null) {
                    all_core_points_map.put(q, new HashMap<>());
                }
                all_core_points_map.get(q).put(sIdx, corepoints.get(q));
            }
        }

        int newestSlide = (currentTime - 1) / Constants.slide;

        for (OD_Query q : all_queries) {
            for (CorePoint c : all_distinct_core_map.get(q)) {
                if (!c.totalHalfRPoints.containsKey(q)) {
                    c.totalHalfRPoints.put(q, c.getTotalHalfRPoints(q));
                } else {
                    c.totalHalfRPoints.put(q, c.totalHalfRPoints.get(q) + c.getTotalHalfRPoints(q, newestSlide));
                }
            }
        }
        
        //sort queries in increasing R order
        Collections.sort(all_queries, (OD_Query o1, OD_Query o2) -> {

            if (o1.R > o2.R) {
                return 2;
            } else if (o1.R == o2.R) {
                if (o1.k < o2.k) {
                    return 1;
                } else if (o1.k > o2.k) {
                    return -1;
                } else {
                    return 0;
                }
            } else {
                return -2;
            }
        });
        
        
        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {

                ArrayList<OD_Query> clone_all_queries = new ArrayList<>();
                for (OD_Query q : all_queries) {
                    clone_all_queries.add(q);
                }
                while (!clone_all_queries.isEmpty()) {
                    ArrayList<OD_Query> can_skip = new ArrayList<>();
//                    int selected_idx = 0;
//                    int selected_idx = clone_all_queries.size() / 2;
                    int selected_idx = clone_all_queries.size() - 1;
//                    Random r = new Random();
//                    int selected_idx = r.nextInt(clone_all_queries.size());
//                    OD_Query q = clone_all_queries.get(clone_all_queries.size() / 2);
                    OD_Query q = clone_all_queries.get(selected_idx);
                    if (d.closeCoreMaps_halfR.get(q) != null && d.closeCoreMaps_halfR.get(q).totalHalfRPoints.get(q) >= q.k + 1) {
                        //inlier
                        for (int t = selected_idx + 1; t < clone_all_queries.size(); t++) {
                            OD_Query q2 = clone_all_queries.get(t);
                            if (q2.k <= q.k) {
                                can_skip.add(q2);
                            }
                        }

                    }
                    if (d.neighborCount.get(q) == null || d.neighborCount.get(q) < q.k) {
                        probe(d, newestSlide, q);
                    }

                    if (d.neighborCount.get(q) < q.k) {
                        //outlier
                        result.get(q).add(d);
                        //will be outlier with smaller r and larger k
                        for (int t = 0; t < selected_idx; t++) {
                            OD_Query q2 = clone_all_queries.get(t);
                            if (q2.k >= q.k) {
                                result.get(q2).add(d);
                                can_skip.add(q2);
                            }
                        }
                    } else {
                        //inlier
                        //will be inlier with larger r and smaller k
                        for (int t = selected_idx + 1; t < clone_all_queries.size(); t++) {
                            OD_Query q2 = clone_all_queries.get(t);
                            if (q2.k <= q.k) {
                                can_skip.add(q2);
                            }
                        }
                    }
                    for (OD_Query q2 : can_skip) {
                        clone_all_queries.remove(q2);

                    }
                    clone_all_queries.remove(q);
                }

            }
        }

        

        return result;
    }

    private ResultFindCore findCloseCore(C_Data d, int slideIndex, OD_Query q) {

        ArrayList<CorePoint> resultCore = null;

        if (d.closeCoreMaps_halfR.get(q) != null && all_core_points_map.get(q).get(slideIndex).contains(d.closeCoreMaps_halfR.get(q)) //                && d.closeCoreMaps_halfR.get(q).closeNeighbors_halfR.containsKey(slideIndex)
                ) {
            resultCore = new ArrayList<>();
            resultCore.add(d.closeCoreMaps_halfR.get(q));
            return new ResultFindCore(q.R / 2, resultCore);
        } else if (d.closeCoreMaps_R.containsKey(q) && !d.closeCoreMaps_R.get(q).isEmpty()) {
            for (CorePoint c : d.closeCoreMaps_R.get(q)) {
                if (all_core_points_map.get(q).get(slideIndex).contains(c)) {
                    resultCore = new ArrayList<>();
                    resultCore.add(c);
                    return new ResultFindCore(q.R, resultCore);
                }
            }

        }

        ArrayList<CorePoint> corePoints = all_core_points_map.get(q).get(slideIndex);

        ArrayList<CorePoint> inRangeRCores = new ArrayList<>();
        ArrayList<CorePoint> inRangeDoubleRCores = new ArrayList<>();
        ArrayList<Double> distance_to_cores = new ArrayList<>();
        if (corePoints != null) {
            for (int i = 0; i < corePoints.size(); i++) {
                CorePoint c = corePoints.get(i);
                double distance = DistanceFunction.euclideanDistance(d, c);

                if (distance <= q.R) {
                    inRangeRCores.add(c);
                    break;
                } else if (distance <= q.R * 2) {
                    inRangeDoubleRCores.add(c);
                    distance_to_cores.add(distance);
                }
            }
        }
        if (!inRangeRCores.isEmpty()) {
            return new ResultFindCore(q.R, inRangeRCores);
        } else if (!inRangeDoubleRCores.isEmpty()) {
            return new ResultFindCore(q.R * 2, inRangeDoubleRCores, distance_to_cores);
        } else {
            return new ResultFindCore(q.R * 2, null);
        }
    }

    private void probe_slide_left(C_Data d, int slideIndex, OD_Query q) {

        int oldNumNeighbor = 0;
        if (d.neighborCount.get(q) != null) {
            oldNumNeighbor = d.neighborCount.get(q);
        }
        //scan possible points 
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();
        //find close core
        ResultFindCore rf = findCloseCore(d, slideIndex, q);

        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        int case_ = 0;
        if (cores != null) {
            if (distance <= q.R / 2) {
                CorePoint c = cores.get(0);
                //grab close neighbor in range R/2 of c
                if (d.neighborCount.get(q) != null) {
                    d.neighborCount.put(q, d.neighborCount.get(q) + c.getTotalHalfRPoints(q, slideIndex));
                } else {
                    d.neighborCount.put(q, c.getTotalHalfRPoints(q, slideIndex));
                }
                if (d.neighborCount.get(q) >= q.k) {

                    if (d.pred_neighbor_count.get(q) == null) {
                        d.pred_neighbor_count.put(q, new HashMap<>());
                    }

                    d.pred_neighbor_count.get(q).put(slideIndex, c.getTotalHalfRPoints(q, slideIndex));

                    if (!neighborCountTrigger.containsKey(q)) {
                        neighborCountTrigger.put(q, new HashMap<>());
                    }
                    if (neighborCountTrigger.get(q).containsKey(slideIndex)) {
                        neighborCountTrigger.get(q).get(slideIndex).add(d);
                    } else {
                        HashSet<C_Data> hs = new HashSet<>();
                        hs.add(d);
                        neighborCountTrigger.get(q).put(slideIndex, hs);
                    }

                    return;

                }

                possibleCandidates.add(c.getDataInRange(q.R / 2, 3 * q.R / 2, slideIndex));

            } else if (distance <= q.R) {

                possibleCandidates.add(cores.get(0).getDataInRange(0, 2 * q.R, slideIndex));

            } else if (distance <= q.R * 2) {
                case_ = 1;
                for (int i = 0; i < cores.size(); i++) {
                    CorePoint c = cores.get(i);
                    if (rf.distance_to_cores.get(i) <= q.R * 3 / 2) {
                        possibleCandidates.add(c.getDataInRange(0, q.R / 2, slideIndex));
                    }

                }

                for (CorePoint c : cores) {
                    possibleCandidates.add(c.getDataInRange(q.R / 2, q.R, slideIndex));

                }

            }

//        start = Utils.getCPUTime();
//        HashSet<C_Data> checked = new HashSet<>();
            int min_arrival_time = all_slides.get(slideIndex).get(0).arrivalTime;

            boolean[] checked = null;
            if (case_ == 1) {
                checked = new boolean[Constants.slide];
            }

            outerloop:
            for (ArrayList<C_Data> ps : possibleCandidates) {

                for (int t = 0; t < ps.size(); t++) {
                    C_Data d2 = ps.get(t);
//                if (!checked.contains(d2)) {
                    if (case_ == 0 || (case_ == 1 && !checked[d2.arrivalTime - min_arrival_time])) {

                        if (DistanceFunction.euclideanDistance(d, d2) <= q.R) {
                            d.neighborCount.put(q, d.neighborCount.get(q) + 1);

                            if (d.neighborCount.get(q) >= q.k) {

                                break outerloop;
                            }
                        }

                        if (case_ == 1) {
                            checked[d2.arrivalTime - min_arrival_time] = true;
                        }
                    }
                }

            }
            if (d.pred_neighbor_count.get(q) == null) {
                d.pred_neighbor_count.put(q, new HashMap<>());
            }
            d.pred_neighbor_count.get(q).put(slideIndex, d.neighborCount.get(q) - oldNumNeighbor);
//            start = Utils.getCPUTime();

            if (neighborCountTrigger.get(q) == null) {
                neighborCountTrigger.put(q, new HashMap<>());
            }
            if (neighborCountTrigger.get(q).containsKey(slideIndex)) {
                neighborCountTrigger.get(q).get(slideIndex).add(d);
            } else {
                HashSet<C_Data> hs = new HashSet<>();
                hs.add(d);
                neighborCountTrigger.get(q).put(slideIndex, hs);
            }
        }

    }

    private void probe_slide_right(C_Data d, int slideIndex, OD_Query q) {

        //scan possible points 
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();

        //find close core
        ResultFindCore rf = findCloseCore(d, slideIndex, q);
//        System.out.println("Time to find core = "+ (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        double distance = rf.getDistance();
        ArrayList<CorePoint> cores = rf.getCore();
        int case_ = 0;
        if (cores != null) {
            if (distance <= q.R / 2) {
                CorePoint c = cores.get(0);

                //grab close neighbor in range R/2 of c
                if (d.neighborCount.get(q) != null) {
                    d.neighborCount.put(q, d.neighborCount.get(q) + c.getTotalHalfRPoints(q, slideIndex));
                } else {
                    d.neighborCount.put(q, c.getTotalHalfRPoints(q, slideIndex));
                }
                if (d.numSucceedingNeighbor.get(q) == null) {
                    d.numSucceedingNeighbor.put(q, c.getTotalHalfRPoints(q, slideIndex));
                } else {
                    d.numSucceedingNeighbor.put(q, d.numSucceedingNeighbor.get(q) + c.getTotalHalfRPoints(q, slideIndex));
                }

                if (d.numSucceedingNeighbor.get(q) >= q.k) {

                    return;

                }

                possibleCandidates.add(c.getDataInRange(q.R / 2, 3 * q.R / 2, slideIndex));

            } else if (distance <= q.R) {

                possibleCandidates.add(cores.get(0).getDataInRange(0, 2 * q.R, slideIndex));

            } else if (distance <= q.R * 2) {
                case_ = 1;
                for (int i = 0; i < cores.size(); i++) {
                    CorePoint c = cores.get(i);
                    if (rf.distance_to_cores.get(i) <= q.R * 3 / 2) {
                        possibleCandidates.add(c.getDataInRange(0, q.R / 2, slideIndex));
                    }

                }
                for (CorePoint c : cores) {
                    possibleCandidates.add(c.getDataInRange(q.R / 2, q.R, slideIndex));
                }

            }

            int min_arrival_time = all_slides.get(slideIndex).get(0).arrivalTime;

            boolean[] checked = null;
            if (case_ == 1) {
                checked = new boolean[Constants.slide];
            }
            int oldNumSucNeighbor = 0;
            if (d.numSucceedingNeighbor.containsKey(q)) {
                oldNumSucNeighbor = d.numSucceedingNeighbor.get(q);
            }
//        outerloop:
            for (ArrayList<C_Data> ps : possibleCandidates) {

                for (int t = 0; t < ps.size(); t++) {
                    C_Data d2 = ps.get(t);
//                if (!checked.contains(d2)) {
                    if (case_ == 0 || (case_ == 1 && !checked[d2.arrivalTime - min_arrival_time])) {

                        if (DistanceFunction.euclideanDistance(d, d2) <= q.R) {
                            if (d.numSucceedingNeighbor.get(q) != null) {
                                d.numSucceedingNeighbor.put(q, d.numSucceedingNeighbor.get(q) + 1);
                            } else {
                                d.numSucceedingNeighbor.put(q, 1);
                            }

                            if (d.numSucceedingNeighbor.get(q) >= q.k) {
                                //test remove preceding neighbor map 
                                if (d.pred_neighbor_count.containsKey(q)) {
                                    d.pred_neighbor_count.get(q).clear();
                                }

                                if (d.neighborCount.containsKey(q)) {
                                    d.neighborCount.put(q, d.neighborCount.get(q) + d.numSucceedingNeighbor.get(q) - oldNumSucNeighbor);
                                } else {
                                    d.neighborCount.put(q, d.numSucceedingNeighbor.get(q) - oldNumSucNeighbor);
                                }
                                return;
                            }
                        }
//                    checked.add(d2);
                        if (case_ == 1) {
                            checked[d2.arrivalTime - min_arrival_time] = true;
                        }
                    }
                }

            }
            if (d.neighborCount.containsKey(q)) {
                d.neighborCount.put(q, d.neighborCount.get(q) + d.numSucceedingNeighbor.get(q) - oldNumSucNeighbor);
            } else {
                d.neighborCount.put(q, d.numSucceedingNeighbor.get(q) - oldNumSucNeighbor);
            }
        } else {
//            System.out.println("No core found!!!");
        }
    }

    private void probe(C_Data d, int newestSlide, OD_Query q) {
        if (d.lastProbRight.get(q) == null || d.lastProbRight.get(q) < newestSlide) {
            //prob right first
            int slideIndex = d.sIndex;
            if (d.lastProbRight.get(q) != null) {
                slideIndex = d.lastProbRight.get(q) + 1;
            }

            while (slideIndex <= newestSlide && (d.neighborCount.get(q) == null || d.neighborCount.get(q) < q.k)) {
                probe_slide_right(d, slideIndex, q);
                d.lastProbRight.put(q, slideIndex);
                slideIndex++;
            }
        }
        //prob left
        if (d.neighborCount.get(q) < q.k) {

            int slideIndex = d.sIndex - 1;
            if (d.lastProbLeft.get(q) != null) {
                slideIndex = d.lastProbLeft.get(q) - 1;
            }

            while (slideIndex > expiredSlideIndex && slideIndex >= 0
                    && (d.neighborCount.get(q) == null || d.neighborCount.get(q) < q.k)) {
                probe_slide_left(d, slideIndex, q);
                d.lastProbLeft.put(q, slideIndex);
                slideIndex--;
            }
        }
//        }
//        long start = Utils.getCPUTime();
        if (d.neighborCount.get(q) < q.k) {
            //add to outlier List
            if (!outlierList_map.containsKey(q)) {
                outlierList_map.put(q, new HashMap<>());
            }
            if (outlierList_map.get(q).containsKey(d.sIndex)) {
                outlierList_map.get(q).get(d.sIndex).add(d);
            } else {
                HashSet hs = new HashSet();
                hs.add(d);
                outlierList_map.get(q).put(d.sIndex, hs);
            }
        }
//        timeForAddingToOutlierList += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
    }

    private HashMap<OD_Query, ArrayList<CorePoint>> selectCoreForAll(Integer sIdx, ArrayList<OD_Query> queries) {
        HashMap<OD_Query, ArrayList<CorePoint>> result = new HashMap<>();
        HashMap<CorePoint, Boolean> probed_core = new HashMap<>();
        for (OD_Query q : queries) {
            result.put(q, new ArrayList<>());
            if (all_mtree_map.get(q) == null) {
                all_mtree_map.put(q, new MTreeCorePoint());
            }
        }
        Collections.sort(queries, (OD_Query o1, OD_Query o2) -> {
            if (o1.R > o2.R) {
                return -1;
            } else if (o1.R == o2.R) {
                return 0;
            } else {
                return 1;
            }
        });
        for (int q_idx = 0; q_idx < queries.size(); q_idx++) {
            OD_Query q = queries.get(q_idx);
            //select core from largest R level
            ArrayList<CorePoint> corePoints = result.get(q);
            ArrayList<CorePoint> newCores = new ArrayList<>();

            for (int i = 0; i < Constants.slide; i++) {

                C_Data d = all_slides.get(sIdx).get(i);

                //check if data is not connected to a core point
                if ((d.closeCoreMaps_R.get(q) == null && d.closeCoreMaps_halfR.get(q) == null)) {
                    //scan with current cores first
                    for (int j = corePoints.size() - 1; j >= 0; j--) {
                        CorePoint c = corePoints.get(j);
                        if (!probed_core.containsKey(c) || probed_core.get(c) == false) {
                            double distance = DistanceFunction.euclideanDistance(d, c);

                            if (distance <= q.R) {
                                for (int t = q_idx; t < queries.size(); t++) {
                                    if (distance <= queries.get(t).R / 2) {
                                        d.closeCoreMaps_halfR.put(queries.get(t), c);
                                    } else if (distance <= queries.get(t).R) {
                                        if (d.closeCoreMaps_R.get(queries.get(t)) == null) {
                                            ArrayList<CorePoint> arr = new ArrayList<>();
                                            d.closeCoreMaps_R.put(queries.get(t), arr);
                                        }
                                        d.closeCoreMaps_R.get(queries.get(t)).add(c);
                                    }

                                }
                                c.putDataToBin(d, distance, sIdx);
                                break;
                            }
                        }

                    }
                }

                if ((d.closeCoreMaps_R.get(q) == null && d.closeCoreMaps_halfR.get(q) == null)) {
                    //using mtree
                    MTreeCorePoint.Query query = all_mtree_map.get(q).getNearest(d, q.R, 1);
                    CorePoint c = null;
                    double distance = Double.MAX_VALUE;
                    for (MTreeClass.ResultItem ri : query) {
                        c = (CorePoint) ri.data;
                        distance = ri.distance;
                    }
                    if (distance <= q.R && c != null) {
//                        needProbCore.add(c);
                        if (!probed_core.containsKey(c) || probed_core.get(c) == false) {
                            c.putDataToBin(d, distance, sIdx);
                        }
                        //add c to the core of this slide and all the queries with smaller R
//                        corePoints.add(c);
                        for (int t = q_idx; t < queries.size(); t++) {
                            result.get(queries.get(t)).add(c);
                            if (distance <= queries.get(t).R) {

                                if (distance <= queries.get(t).R / 2) {

                                    d.closeCoreMaps_halfR.put(queries.get(t), c);

                                } else {
                                    if (d.closeCoreMaps_R.get(queries.get(t)) == null) {
                                        ArrayList<CorePoint> arr = new ArrayList<>();
                                        d.closeCoreMaps_R.put(queries.get(t), arr);

                                    }
                                    d.closeCoreMaps_R.get(queries.get(t)).add(c);

                                }
                            }
                        }

//                    scanForCore(c, sIdx);
                    } else {
                        c = new CorePoint(d);
//                        needProbCore.add(c);
                        newCores.add(c);

                        for (int t = q_idx; t < queries.size(); t++) {
                            result.get(queries.get(t)).add(c);

                            d.closeCoreMaps_halfR.put(queries.get(t), c);

                        }
                        c.creatBins(queries.subList(q_idx, queries.size()));
                        c.putDataToBin(d, 0.0, sIdx);

                    }

                }
            }

            //find scan for cores
            boolean[] checked = new boolean[Constants.slide];

            for (CorePoint c : corePoints) {
                if (!probed_core.containsKey(c) || probed_core.get(c) == false) {
                    for (int i = 0; i < Constants.slide; i++) {
                        checked[i] = false;
                    }
                    for (CorePoint c2 : corePoints) {

                        if (c != c2) {
                            double distance = DistanceFunction.euclideanDistance(c, c2);

                            if (distance <= q.R * 3) {
                                checked = probCoreWithList(c, c2.getDataInRange(0, q.R, sIdx), sIdx, queries.subList(q_idx, queries.size()),
                                        checked, all_slides.get(sIdx).get(0).arrivalTime);
                            }
                        }
                    }
                    probed_core.put(c, true);
                }

            }
            for (CorePoint c : newCores) {
                all_distince_core_point.add(c);
                for (int t = q_idx; t < queries.size(); t++) {
                    all_mtree_map.get(queries.get(t)).add(c);
                    if (!all_distinct_core_map.containsKey(queries.get(t))) {
                        all_distinct_core_map.put(queries.get(t), new ArrayList<>());
                    }
                    all_distinct_core_map.get(queries.get(t)).add(c);
                }

            }
        }

        return result;

    }

    private boolean[] probCoreWithList(CorePoint c, ArrayList<C_Data> candidates, int sIdx, List<OD_Query> queries,
            boolean[] checked, int start_time) {

        //create bins if not exists
        c.creatBins(queries);

        if (candidates != null) {
            for (C_Data d2 : candidates) {
                if (!checked[d2.arrivalTime - start_time]) {
                    double distance = DistanceFunction.euclideanDistance(c, d2);
                    if (distance <= queries.get(0).R * 2) //put d2 to correct bin
                    {
                        c.putDataToBin(d2, distance, sIdx);

                        for (OD_Query q : queries) {
                            if (distance <= q.R / 2) {
                                d2.closeCoreMaps_halfR.put(q, c);
                            } else if (distance <= q.R) {
                                if (d2.closeCoreMaps_R.get(q) == null) {
                                    d2.closeCoreMaps_R.put(q, new ArrayList<>());
                                }
                                d2.closeCoreMaps_R.get(q).add(c);
                            }
                        }
                    }
                    checked[d2.arrivalTime - start_time] = true;
                }
            }
        }
        return checked;
    }

    private void processExpiredData(int expiredSlideIndex) {
        all_slides.remove(expiredSlideIndex);
        for (OD_Query q : all_queries) {
            HashMap<Integer, HashSet<C_Data>> outlierList = outlierList_map.get(q);
            if (outlierList != null) {
                outlierList.remove(expiredSlideIndex);
            }

            if (neighborCountTrigger.get(q) != null) {
                if (neighborCountTrigger.get(q).containsKey(expiredSlideIndex)) {
                    for (C_Data d : neighborCountTrigger.get(q).get(expiredSlideIndex)) {
                        if (d.pred_neighbor_count.get(q).containsKey(expiredSlideIndex)) {
                            d.neighborCount.put(q, d.neighborCount.get(q) - d.pred_neighbor_count.get(q).get(expiredSlideIndex));
                            d.pred_neighbor_count.get(q).remove(expiredSlideIndex);
                        }
                    }
                }
                if (neighborCountTrigger.get(q).containsKey(expiredSlideIndex - 1)) {
                    neighborCountTrigger.get(q).remove(expiredSlideIndex - 1);
                }
            }
            if (all_core_points_map.get(q) != null) {
                all_core_points_map.get(q).remove(expiredSlideIndex);
            }
        }

        for (CorePoint c : all_distince_core_point) {
            for (Bin b : c.all_bins) {
                if (b.data.containsKey(expiredSlideIndex)) {
                    for (OD_Query q : all_queries) {
                        if (c.totalHalfRPoints.containsKey(q)) {
                            if (b.max_val <= q.R / 2) {
                                c.totalHalfRPoints.put(q, c.totalHalfRPoints.get(q) - b.data.get(expiredSlideIndex).size());
                            }
                        }
                    }
                    b.data.remove(expiredSlideIndex);

                }
            }

        }

    }

    class Bin {

        public double min_val;
        public double max_val;
        public HashMap<Integer, ArrayList<C_Data>> data;

        public Bin(double min, double max) {
            this.min_val = min;
            this.max_val = max;
            this.data = new HashMap<>();
        }
    }

    class CorePoint extends C_Data implements Comparable<Data> {

        public ArrayList<Bin> all_bins = new ArrayList<>();

        public ArrayList<C_Data> getDataInRange(double min, double max) {
            ArrayList<C_Data> result = new ArrayList<>();
            for (Bin b : all_bins) {
                if (b.min_val >= min && b.max_val <= max) {
                    for (Integer sIdx : b.data.keySet()) {
                        result.addAll(b.data.get(sIdx));
                    }
                }
            }
            return result;
        }

        public void creatBins(List<OD_Query> queries) {
            //create bins if not exists
            if (all_bins == null || all_bins.isEmpty()) {
                all_bins = new ArrayList<>();
                ArrayList<Double> all_r = new ArrayList<>();
                all_r.add(new Double(0));
                for (OD_Query q : queries) {
                    all_r.add(q.R / 2);
                    all_r.add(q.R);
                    all_r.add(3 * q.R / 2);
                    all_r.add(2 * q.R);
                }
                Collections.sort(all_r);
                for (int i = 1; i < all_r.size(); i++) {
                    Bin b = new Bin(all_r.get(i - 1), all_r.get(i));
                    all_bins.add(b);
                }
            }
        }

        public void putDataToBin(C_Data d, double distance, int sIdx) {
            for (Bin b : all_bins) {
                if ((b.min_val <= distance) && b.max_val >= distance) {
                    if (!b.data.containsKey(sIdx)) {
                        b.data.put(sIdx, new ArrayList<>());
                    }
                    b.data.get(sIdx).add(d);
                }
            }
        }

        public ArrayList<C_Data> getDataInRange(double min, double max, int sIdx) {
            ArrayList<C_Data> result = new ArrayList<>();
            for (Bin b : all_bins) {
                if (b.min_val >= min && b.max_val <= max) {
                    if (b.data.containsKey(sIdx)) {
                        result.addAll(b.data.get(sIdx));
                    }
                }
            }
            return result;
        }

        public HashMap<OD_Query, Integer> totalHalfRPoints = new HashMap<>();

        public int getTotalHalfRPoints(OD_Query q) {
            int t = 0;
            for (Bin b : all_bins) {
                if (b.max_val <= q.R / 2) {
                    for (Integer idx : b.data.keySet()) {

                        t += b.data.get(idx).size();

                    }
                }
            }
            return t;
        }

        public int getTotalHalfRPoints(OD_Query q, int sIdx) {
            int t = 0;
            for (Bin b : all_bins) {
                if (b.max_val <= q.R / 2) {
                    if (b.data.containsKey(sIdx)) {
                        t += b.data.get(sIdx).size();
                    }

                }
            }
            return t;
        }

//        public int getTotal32RPoints(OD_Query q) {
//            int t = 0;
//            for (Map.Entry<Integer, ArrayList<C_Data>> e : closeNeighbors_3halfR.get(q).entrySet()) {
//                t += e.getValue().size();
//            }
//            return t;
//        }
//
//        public int getTotal2RPoints(OD_Query q) {
//            int t = 0;
//            for (Map.Entry<Integer, ArrayList<C_Data>> e : closeNeighbors_2R.get(q).entrySet()) {
//                t += e.getValue().size();
//            }
//            return t;
//        }
//
//        public int getTotalRPoints(OD_Query q) {
//            int t = 0;
//            for (Map.Entry<Integer, ArrayList<C_Data>> e : closeNeighbors_2R.get(q).entrySet()) {
//                t += e.getValue().size();
//            }
//            return t;
//        }
//
//        public boolean isCoveredAllSlides(OD_Query q) {
//            for (Map.Entry<Integer, ArrayList<C_Data>> e : all_slides.entrySet()) {
//                if (!closeNeighbors_halfR.get(q).containsKey(e.getKey())) {
//                    return false;
//                }
//            }
//            return true;
//        }
        public CorePoint(C_Data d) {
            this.values = d.values;
            this.hashCode = d.hashCode;
            this.arrivalTime = d.arrivalTime;

        }

    }

    class C_Data extends Data {

        private HashMap<OD_Query, Integer> numSucceedingNeighbor = new HashMap<>();

//        private boolean isOutlier;
        public HashMap<OD_Query, Integer> lastProbRight = new HashMap<>();
        public HashMap<OD_Query, Integer> lastProbLeft = new HashMap<>();

        private HashMap<OD_Query, HashMap<Integer, Integer>> pred_neighbor_count = new HashMap<>();
        public HashMap<OD_Query, Integer> neighborCount = new HashMap<>();
//
        public HashMap<OD_Query, CorePoint> closeCoreMaps_halfR = new HashMap<>();
        public HashMap<OD_Query, ArrayList<CorePoint>> closeCoreMaps_R = new HashMap<>();

        public int sIndex = -1;

        public HashMap<OD_Query, Boolean> status = new HashMap<>();

        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;

            this.sIndex = (arrivalTime - 1) / Constants.slide;
        }

        public C_Data() {

        }

//        public int countNeighbor() {
//            return neighborCount;
//        }
    }

//    public class OD_Query {
//
//        public double R;
//        public int k;
////        public int id;
//
//        public OD_Query(double R_, int k_) {
//            this.R = R_;
//            this.k = k_;
//        }
//    }
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

}
