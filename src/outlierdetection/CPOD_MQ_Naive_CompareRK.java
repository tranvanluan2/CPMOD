/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import mtree.tests.Data;
import mtree.utils.Constants;

/**
 *
 * @author luan
 */
public class CPOD_MQ_Naive_CompareRK {

    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public static HashMap<OD_Query, HashMap<Integer, ArrayList<CorePoint>>> all_core_points_map = new HashMap<>();
    public static HashMap<OD_Query, MTreeCorePoint> all_mtree_map = new HashMap<>();

    public static HashMap<OD_Query, ArrayList<CorePoint>> all_distinct_core_map = new HashMap<>();
//    public static HashMap<OD_Query, HashMap<Integer, HashSet<C_Data>>> outlierList_map = new HashMap<>();

    public ArrayList<OD_Query> all_queries = new ArrayList<>();

//    public static HashMap<OD_Query, HashMap<Integer, HashSet<C_Data>>> neighborCountTrigger = new HashMap<>();
    public long numDCS = 0;

    public void add_query(OD_Query q) {
        all_queries.add(q);

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
        int newestSlide = (currentTime - 1) / Constants.slide;
        processExpiredData(expiredSlideIndex, newestSlide);

        for (OD_Query q : all_queries) {
            int query_expire_slide = newestSlide - q.W / Constants.slide;
            for (int sIdx : slide_to_process) {
                if (sIdx > query_expire_slide) {
                    ArrayList<CorePoint> corePoints = selectCore(sIdx, q);
                    if (all_core_points_map.get(q) == null) {
                        all_core_points_map.put(q, new HashMap<>());
                    }
                    all_core_points_map.get(q).put(sIdx, corePoints);
                }
            }
        }

        if (currentTime == Constants.W) {
            for (OD_Query q : all_queries) {
                for (CorePoint c : all_distinct_core_map.get(q)) {
                    c.totalHalfRPoints = c.getTotalHalfRPoints();
                }
            }
        } else if (data.size() == Constants.slide) {
            for (OD_Query q : all_queries) {
                for (CorePoint c : all_distinct_core_map.get(q)) {
                    if (c.closeNeighbors_halfR.get(newestSlide) != null) {
                        c.totalHalfRPoints += c.closeNeighbors_halfR.get(newestSlide).size();
                    }
                }
            }
        }

        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {

                ArrayList<OD_Query> clone_all_queries = new ArrayList<>();
                for (OD_Query q : all_queries) {

                    clone_all_queries.add(q);
                }

                while (!clone_all_queries.isEmpty()) {
                    int selected_idx = clone_all_queries.size() - 1;
                    OD_Query q = clone_all_queries.get(selected_idx);
                    if (currentTime - q.W < d.arrivalTime) {
                        if (d.closeCoreMaps_halfR.get(q) != null
                                && d.closeCoreMaps_halfR.get(q).totalHalfRPoints >= q.k + 1) {

                        } else if (d.neighborCount.get(q) == null || d.neighborCount.get(q) < q.k) {
                            probe(d, newestSlide, q);
                            if (d.neighborCount.get(q) == null || d.neighborCount.get(q) < q.k) {
                                if ((currentTime - Constants.W) / Constants.slide * Constants.slide == (currentTime - Constants.W)) //outlier
                                {
                                    result.get(q).add(d);
                                }
                                //will be outlier with smaller r and larger k
                            }
                        }
                    }

                    clone_all_queries.remove(q);
                }

            }
        }

        return result;
    }

    private ResultFindCore findCloseCore(C_Data d, int slideIndex, OD_Query q) {

        ArrayList<CorePoint> resultCore = null;

        if (d.closeCoreMaps_halfR.get(q) != null
                && d.closeCoreMaps_halfR.get(q).closeNeighbors_halfR.containsKey(slideIndex)) {
            resultCore = new ArrayList<>();
            resultCore.add(d.closeCoreMaps_halfR.get(q));
            return new ResultFindCore(q.R / 2, resultCore);
        } else if (d.closeCoreMaps_R.containsKey(q) && !d.closeCoreMaps_R.get(q).isEmpty()) {
            for (CorePoint c : d.closeCoreMaps_R.get(q)) {
                if (c.closeNeighbors_2R.containsKey(slideIndex)) {
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
                    d.neighborCount.put(q, d.neighborCount.get(q) + c.closeNeighbors_halfR.get(slideIndex).size());
                } else {
                    d.neighborCount.put(q, c.closeNeighbors_halfR.get(slideIndex).size());
                }
                if (d.neighborCount.get(q) >= q.k) {

                    if (d.pred_neighbor_count.get(q) == null) {
                        d.pred_neighbor_count.put(q, new HashMap<>());
                    }

                    d.pred_neighbor_count.get(q).put(slideIndex, c.closeNeighbors_halfR.get(slideIndex).size());

//                    if (!neighborCountTrigger.containsKey(q)) {
//                        neighborCountTrigger.put(q, new HashMap<>());
//                    }
//                    if (neighborCountTrigger.get(q).containsKey(slideIndex)) {
//                        neighborCountTrigger.get(q).get(slideIndex).add(d);
//                    } else {
//                        HashSet<C_Data> hs = new HashSet<>();
//                        hs.add(d);
//                        neighborCountTrigger.get(q).put(slideIndex, hs);
//                    }
                    return;

                }
                possibleCandidates.add(c.closeNeighbors_R.get(slideIndex));
                possibleCandidates.add(c.closeNeighbors_3halfR.get(slideIndex));

            } else if (distance <= q.R) {

                possibleCandidates.add(cores.get(0).closeNeighbors_halfR.get(slideIndex));
                possibleCandidates.add(cores.get(0).closeNeighbors_R.get(slideIndex));
                possibleCandidates.add(cores.get(0).closeNeighbors_3halfR.get(slideIndex));
                possibleCandidates.add(cores.get(0).closeNeighbors_2R.get(slideIndex));

            } else if (distance <= q.R * 2) {
                case_ = 1;
                for (int i = 0; i < cores.size(); i++) {
                    CorePoint c = cores.get(i);
                    if (rf.distance_to_cores.get(i) <= q.R * 3 / 2) {
                        possibleCandidates.add(c.closeNeighbors_halfR.get(slideIndex));
                    }

                }

                for (CorePoint c : cores) {
                    possibleCandidates.add(c.closeNeighbors_R.get(slideIndex));

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

                        double dist;
                        if (d2.last_checked_dist != null && d2.last_checked_dist.arrivalTime == d.arrivalTime) {
                            dist = d2.last_distance;
                        } else {
                            dist = DistanceFunction.euclideanDistance(d, d2);
                            d2.last_checked_dist = d;
                            d2.last_distance = dist;
                            numDCS += 1;
                        }

//                        double dist = DistanceFunction.euclideanDistance(d, d2);
//                        numDCS +=1;
                        if (dist <= q.R) {
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

//            if (neighborCountTrigger.get(q) == null) {
//                neighborCountTrigger.put(q, new HashMap<>());
//            }
//            if (neighborCountTrigger.get(q).containsKey(slideIndex)) {
//                neighborCountTrigger.get(q).get(slideIndex).add(d);
//            } else {
//                HashSet<C_Data> hs = new HashSet<>();
//                hs.add(d);
//                neighborCountTrigger.get(q).put(slideIndex, hs);
//            }
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
                    d.neighborCount.put(q, d.neighborCount.get(q) + c.closeNeighbors_halfR.get(slideIndex).size());
                } else {
                    d.neighborCount.put(q, c.closeNeighbors_halfR.get(slideIndex).size());
                }
                if (d.numSucceedingNeighbor.get(q) == null) {
                    d.numSucceedingNeighbor.put(q, c.closeNeighbors_halfR.get(slideIndex).size());
                } else {
                    d.numSucceedingNeighbor.put(q, d.numSucceedingNeighbor.get(q) + c.closeNeighbors_halfR.get(slideIndex).size());
                }

                if (d.numSucceedingNeighbor.get(q) >= q.k) {

                    return;

                }
                possibleCandidates.add(c.closeNeighbors_R.get(slideIndex));
                possibleCandidates.add(c.closeNeighbors_3halfR.get(slideIndex));

            } else if (distance <= q.R) {

                possibleCandidates.add(cores.get(0).closeNeighbors_halfR.get(slideIndex));
                possibleCandidates.add(cores.get(0).closeNeighbors_R.get(slideIndex));
                possibleCandidates.add(cores.get(0).closeNeighbors_3halfR.get(slideIndex));
                possibleCandidates.add(cores.get(0).closeNeighbors_2R.get(slideIndex));

            } else if (distance <= q.R * 2) {
                case_ = 1;
                for (int i = 0; i < cores.size(); i++) {
                    CorePoint c = cores.get(i);
                    if (rf.distance_to_cores.get(i) <= q.R * 3 / 2) {
                        possibleCandidates.add(c.closeNeighbors_halfR.get(slideIndex));
                    }

                }
                for (CorePoint c : cores) {
                    possibleCandidates.add(c.closeNeighbors_R.get(slideIndex));
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
            } else {
                d.numSucceedingNeighbor.put(q, 0);
            }
//        outerloop:
            for (ArrayList<C_Data> ps : possibleCandidates) {

                for (int t = 0; t < ps.size(); t++) {
                    C_Data d2 = ps.get(t);
//                if (!checked.contains(d2)) {
                    if ((case_ == 0 || (case_ == 1 && !checked[d2.arrivalTime - min_arrival_time]))) {

                        double dist;
                        if (d2.last_checked_dist != null && d2.last_checked_dist.arrivalTime == d.arrivalTime) {
                            dist = d2.last_distance;
                        } else {
                            dist = DistanceFunction.euclideanDistance(d, d2);
                            d2.last_checked_dist = d;
                            d2.last_distance = dist;
                            numDCS += 1;
                        }

//                        double dist = DistanceFunction.euclideanDistance(d, d2);
//                        numDCS +=1;
                        if (dist <= q.R) {
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
        if (d.neighborCount.get(q) == null || d.neighborCount.get(q) < q.k) {

            int slideIndex = d.sIndex - 1;
            if (d.lastProbLeft.get(q) != null) {
                slideIndex = d.lastProbLeft.get(q) - 1;
            }

            while (slideIndex > newestSlide - q.W / q.S && slideIndex >= 0
                    && (d.neighborCount.get(q) == null || d.neighborCount.get(q) < q.k)) {
                probe_slide_left(d, slideIndex, q);
                d.lastProbLeft.put(q, slideIndex);
                slideIndex--;
            }
        }
    }

    private ArrayList<CorePoint> selectCore(Integer sIdx, OD_Query q) {

        ArrayList<CorePoint> corePoints = new ArrayList<>();
        ArrayList<CorePoint> newCores = new ArrayList<>();

        for (int i = 0; i < Constants.slide; i++) {

            C_Data d = all_slides.get(sIdx).get(i);
            //scan with current cores first
            for (int j = corePoints.size() - 1; j >= 0; j--) {
                CorePoint c = corePoints.get(j);
                double distance = DistanceFunction.euclideanDistance(d, c);

                if (distance <= q.R / 2) {
                    ArrayList<C_Data> closeNeighbors = c.closeNeighbors_halfR.get(sIdx);
                    if (closeNeighbors == null) {
                        closeNeighbors = new ArrayList<>(Arrays.asList(d));
                        c.closeNeighbors_halfR.put(sIdx, closeNeighbors);
                    } else {
                        closeNeighbors.add(d);
                    }
                    d.closeCoreMaps_halfR.put(q, c);
                    break;
                } else if (distance <= q.R) {
                    ArrayList<C_Data> closeNeighbors = c.closeNeighbors_R.get(sIdx);
                    if (closeNeighbors == null) {
                        closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                        closeNeighbors.add(d);
                        c.closeNeighbors_R.put(sIdx, closeNeighbors);
                    } else {
                        closeNeighbors.add(d);
                    }
                    if (d.closeCoreMaps_R.get(q) == null) {
                        ArrayList<CorePoint> arr = new ArrayList<>();
                        d.closeCoreMaps_R.put(q, arr);

                    }
                    d.closeCoreMaps_R.get(q).add(c);
                    break;
                }
            }

            if (all_mtree_map.get(q) == null) {
                all_mtree_map.put(q, new MTreeCorePoint());
            }

//            System.out.println(d.closeCoreMaps_R.get(q));
//            System.out.println(d.closeCoreMaps_halfR.get(q));
            if ((d.closeCoreMaps_R.get(q) == null && d.closeCoreMaps_halfR.get(q) == null)) {
                //using mtree
                MTreeCorePoint.Query query = all_mtree_map.get(q).getNearest(d, q.R, 1);
                CorePoint c = null;
                double distance = Double.MAX_VALUE;
                for (MTreeClass.ResultItem ri : query) {
                    c = (CorePoint) ri.data;
                    distance = ri.distance;
                }
                if (distance <= q.R) {
                    //add c to the core of this slide
                    corePoints.add(c);

                    if (distance <= q.R / 2) {
                        ArrayList<C_Data> closeNeighbors = c.closeNeighbors_halfR.get(sIdx);
                        if (closeNeighbors == null) {
                            closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                            closeNeighbors.add(d);
                            c.closeNeighbors_halfR.put(sIdx, closeNeighbors);
                        } else {
                            closeNeighbors.add(d);
                        }
                        d.closeCoreMaps_halfR.put(q, c);

                    } else {
                        ArrayList<C_Data> closeNeighbors = c.closeNeighbors_R.get(sIdx);
                        if (closeNeighbors == null) {
                            closeNeighbors = new ArrayList<>(Arrays.asList(d));
//                            closeNeighbors.add(d);
                            c.closeNeighbors_R.put(sIdx, closeNeighbors);
                        } else {
                            closeNeighbors.add(d);
                        }
                        if (d.closeCoreMaps_R.get(q) == null) {
                            ArrayList<CorePoint> arr = new ArrayList<>();
                            d.closeCoreMaps_R.put(q, arr);

                        }
                        d.closeCoreMaps_R.get(q).add(c);

                    }
//                    scanForCore(c, sIdx);

                } else {
                    c = new CorePoint(d);
                    if (all_distinct_core_map.get(q) == null) {
                        all_distinct_core_map.put(q, new ArrayList<>());
                    }
                    all_distinct_core_map.get(q).add(c);
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
                    d.closeCoreMaps_halfR.put(q, c);

                    //probe neighbors for c
//                    scanForCore(c, sIdx);
                }

            }
        }
        boolean[] checked = new boolean[Constants.slide];
        //find scan for cores
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

            for (int i = 0; i < Constants.slide; i++) {
                checked[i] = false;
            }
            for (C_Data d : c.closeNeighbors_halfR.get(sIdx)) {
                checked[d.arrivalTime - all_slides.get(sIdx).get(0).arrivalTime] = true;
            }
            for (C_Data d : c.closeNeighbors_R.get(sIdx)) {
                checked[d.arrivalTime - all_slides.get(sIdx).get(0).arrivalTime] = true;
            }
            for (CorePoint c2 : corePoints) {

                if (c != c2) {
                    double distance = DistanceFunction.euclideanDistance(c, c2);

                    if (distance <= q.R * 3) {
                        checked = probCoreWithList(c, c2.closeNeighbors_halfR.get(sIdx), sIdx, q, checked, all_slides.get(sIdx).get(0).arrivalTime);
                        checked = probCoreWithList(c, c2.closeNeighbors_R.get(sIdx), sIdx, q, checked, all_slides.get(sIdx).get(0).arrivalTime);
                    }
                }
            }

        }
        for (CorePoint c : newCores) {
            all_mtree_map.get(q).add(c);
        }
        return corePoints;

    }

    private boolean[] probCoreWithList(CorePoint c, ArrayList<C_Data> candidates, int sIdx, OD_Query q, boolean[] checked, int start_time) {
        if (candidates != null) {
            for (C_Data d2 : candidates) {
                if (!checked[d2.arrivalTime - start_time]) {
                    double distance = DistanceFunction.euclideanDistance(c, d2);

                    if (distance <= q.R / 2) {
                        c.closeNeighbors_halfR.get(sIdx).add(d2);
                        d2.closeCoreMaps_halfR.put(q, c);
                    } else if (distance <= q.R) {
                        c.closeNeighbors_R.get(sIdx).add(d2);
                        if (d2.closeCoreMaps_R.get(q) == null) {
                            d2.closeCoreMaps_R.put(q, new ArrayList<>());
                        }
                        d2.closeCoreMaps_R.get(q).add(c);
                    } else if (distance <= q.R * 1.5) {
                        c.closeNeighbors_3halfR.get(sIdx).add(d2);

                    } else if (distance <= q.R * 2) {
                        c.closeNeighbors_2R.get(sIdx).add(d2);

                    }
                    checked[d2.arrivalTime - start_time] = true;
                }

            }
//        }

//            c.closeNeighbors_halfR.put(sIdx, neighborsInHalfR);
//            c.closeNeighbors_R.put(sIdx, neighborsInR);
//            c.closeNeighbors_3halfR.put(sIdx, neighborsIn3HalfR);
//            c.closeNeighbors_2R.put(sIdx, neighborsIn2R);
        }
        return checked;
    }

    private void processExpiredData(int expiredSlideIndex, int newestSlide) {
        all_slides.remove(expiredSlideIndex);

        for (Integer sIdx : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIdx)) {
                for (OD_Query q : all_queries) {
                    int query_expire_slide = newestSlide - q.W / Constants.slide;
                    if (d.pred_neighbor_count.containsKey(q)
                            && d.pred_neighbor_count.get(q).containsKey(query_expire_slide)) {
                        d.neighborCount.put(q, d.neighborCount.get(q)
                                - d.pred_neighbor_count.get(q).get(query_expire_slide));
                        d.pred_neighbor_count.get(q).remove(query_expire_slide);
                    }

                }

            }
        }
        for (OD_Query q : all_queries) {
            int query_expire_slide = newestSlide - q.W / Constants.slide;
            if (all_core_points_map.get(q) != null) {
                all_core_points_map.get(q).remove(query_expire_slide);
            }

            if (all_distinct_core_map.get(q) != null) {
                for (CorePoint c : all_distinct_core_map.get(q)) {
                    if (c.closeNeighbors_halfR.get(query_expire_slide) != null) {
                        c.totalHalfRPoints -= c.closeNeighbors_halfR.get(query_expire_slide).size();
                        c.closeNeighbors_halfR.remove(query_expire_slide);
                    }
                    c.closeNeighbors_R.remove(query_expire_slide);
                    c.closeNeighbors_3halfR.remove(query_expire_slide);
                    c.closeNeighbors_2R.remove(query_expire_slide);
                }
            }
        }
    }

//    public static HashMap<Integer, HashSet<C_Data>> neighborCountTrigger = new HashMap<>();
    class CorePoint extends C_Data implements Comparable<Data> {

        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_halfR = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_R = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_3halfR = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_2R = new HashMap<>();

        public int totalHalfRPoints = 0;

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

        public C_Data last_checked_dist = null;
        public double last_distance = -1;

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
