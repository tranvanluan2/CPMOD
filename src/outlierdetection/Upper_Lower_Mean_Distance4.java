/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class Upper_Lower_Mean_Distance4 {
    
    public static int currentTime;
    
    public static int expiredSlideIndex = -1;
    
    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();
    
    public static ArrayList<CorePoint> all_core_points = new ArrayList<>();
    
    public static HashMap<Integer, HashSet<C_Data>> outlierList = new HashMap<>();
    
    public static HashMap<Integer, HashSet<CorePoint>> corePointMapTrigger = new HashMap<>();
    public static HashMap<Integer, HashSet<C_Data>> neighborCountTrigger = new HashMap<>();
    
    public double timeFindingCore = 0;
    public double timeCheckingCandidates = 0;
    public double timeAddingToNeighborCount = 0;
    public double timeMarkProbing = 0;
    public double avg_candidate_size = 0;
    
    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {
        ArrayList<Data> result = new ArrayList<>((int) (data.size() * 1.0 / 100));
        currentTime = _currentTime;
        expiredSlideIndex = (currentTime - 1) / slide - Constants.W / slide;
        System.out.println("Expire slide index = " + expiredSlideIndex);
        
        processExpiredData(expiredSlideIndex);
        
        ArrayList<C_Data> d_to_process = new ArrayList<>(data.size());
//        HashSet<Integer> slide_to_process = new HashSet<>();

        int[] slide_to_process;
        if (data.size() == W) {
            slide_to_process = new int[W / slide];
            for (int i = 0; i < slide_to_process.length; i++) {
                slide_to_process[i] = i;
            }
        } else {
            slide_to_process = new int[]{(currentTime - 1) / slide};
        }
        
        for (int i = 0; i < data.size(); i++) {
            Data o = data.get(i);
            C_Data d = new C_Data(o);
//            System.out.println("Slide index of d = "+ d.sIndex);
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

//        selectCore();
//
//        System.out.println("Num new core = " + all_core_points.size());
        for (Integer sIndex : slide_to_process) {
            matchCore(sIndex);
        }
        System.out.println("Num core after matching core = " + all_core_points.size());
        
        HashMap<CorePoint, ArrayList<C_Data>> core_cdata_maps = notifyCloseCore();
//        HashMap<CorePoint, ArrayList<C_Data>> core_cdata_maps = notify_inR_Core();
        
        System.out.println("Num cores notified = " + core_cdata_maps.size());
        
        
        scanForCore(core_cdata_maps.keySet());

//        ArrayList<C_Data> points = findPointNeedProbe(HashMap<CorePoint, ArrayList<C_Data>>);
//
//        System.out.println("After scanning, num data points need probe = " + points.size());
        System.out.println("Probing");
        
        int newestSlide = (currentTime - 1) / Constants.slide;
        
        for (CorePoint c : core_cdata_maps.keySet()) {
            for (C_Data d : core_cdata_maps.get(c)) {
                probe(d, c, newestSlide);
                if (d.neighborCount < Constants.k) {
                    result.add(d);
                }
            }
        }
        
        for (CorePoint c : all_core_points) {
            if (c.arrivalTime > expiredSlideIndex) {
                if (c.neighborCount < Constants.k) {
                    result.add(c);
                }
            }
        }
        
        return result;
        
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
        
        for (CorePoint c : all_core_points) {
            c.closeNeighbors_halfR.remove(expiredSlideIndex);
            c.closeNeighbors_R.remove(expiredSlideIndex);
            c.closeNeighbors_3halfR.remove(expiredSlideIndex);
            c.closeNeighbors_2R.remove(expiredSlideIndex);
        }

//        all_core_points.remove(expiredSlideIndex);
    }
    
    private void findCloseCore(C_Data d) {
        int sIndex = d.sIndex;
        //prob for finding close core
        for (CorePoint c : all_core_points) {
            double distance = DistanceFunction.euclideanDistance(c, d);
            
            if (distance <= Constants.R / 2) {
                d.closeCoreMaps_halfR.add(c);
                c.neighborCount += 1;
                if (c.closeNeighbors_halfR.containsKey(sIndex)) {
                    c.closeNeighbors_halfR.get(sIndex).add(d);
                } else {
                    ArrayList<C_Data> neighbors = new ArrayList<>();
                    neighbors.add(d);
                    c.closeNeighbors_halfR.put(sIndex, neighbors);
                }
                break;
                
            } else if (distance <= Constants.R) {
                d.closeCoreMaps_R.add(c);
                c.neighborCount += 1;
                if (c.closeNeighbors_R.containsKey(sIndex)) {
                    c.closeNeighbors_R.get(sIndex).add(d);
                } else {
                    ArrayList<C_Data> neighbors = new ArrayList<>();
                    neighbors.add(d);
                    c.closeNeighbors_R.put(sIndex, neighbors);
                }
//                break;
            }
//            else if (distance <= Constants.R * 1.5) {
//                c.neighborCount += 1;
//                if (c.closeNeighbors_3halfR.containsKey(sIndex)) {
//                    c.closeNeighbors_3halfR.get(sIndex).add(d);
//                } else {
//                    ArrayList<C_Data> neighbors = new ArrayList<>();
//                    neighbors.add(d);
//                    c.closeNeighbors_3halfR.put(sIndex, neighbors);
//                }
//            } else if (distance <= Constants.R * 2) {
//                c.neighborCount += 1;
//                if (c.closeNeighbors_2R.containsKey(sIndex)) {
//                    c.closeNeighbors_2R.get(sIndex).add(d);
//                } else {
//                    ArrayList<C_Data> neighbors = new ArrayList<>();
//                    neighbors.add(d);
//                    c.closeNeighbors_2R.put(sIndex, neighbors);
//                }
//            }
        }
    }
    
    private void selectCore() {
        for (Integer sIndex : all_slides.keySet()) {
            
            matchCore(sIndex);
            
        }
    }
    
    private void matchCore(int sIndex) {
        for (int i = 0; i < Constants.slide; i++) {
            C_Data d = all_slides.get(sIndex).get(i);
            findCloseCore(d);
            if (d.closeCoreMaps_R.isEmpty() && d.closeCoreMaps_halfR.isEmpty()) {
                //promote d to be core
                CorePoint c = new CorePoint(d);
                
                all_slides.get(sIndex).set(i, c);
                all_core_points.add(c);
                
            }
        }
    }
    
    private HashMap<CorePoint, ArrayList<C_Data>> notifyCore() {
        HashMap<CorePoint, ArrayList<C_Data>> result = new HashMap<>();
        for (CorePoint c : all_core_points) {
            if (c.neighborCount < Constants.k) {
                result.put(c, new ArrayList<>());
            }
        }
        for (Integer sIndex : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIndex)) {
                if (!(d instanceof CorePoint)) {
                    if (d.neighborCount < Constants.k) {
                        boolean needCheck = true;
                        if (!d.closeCoreMaps_halfR.isEmpty()) {
                            for (CorePoint c : d.closeCoreMaps_halfR) {
                                if (c.totalCloseNeighbor() >= Constants.k) {
                                    needCheck = false;
                                }
                            }
                        }
                        
                        if (needCheck) {
                            boolean closeCoreChecked = false;
                            for (CorePoint c : d.closeCoreMaps_halfR) {
                                if (result.containsKey(c)) {
                                    closeCoreChecked = true;
                                    result.get(c).add(d);
                                    break;
                                }
                                
                            }
                            if (!closeCoreChecked) {
                                for (CorePoint c : d.closeCoreMaps_R) {
                                    if (result.containsKey(c)) {
                                        closeCoreChecked = true;
                                        result.get(c).add(d);
                                        break;
                                    }
                                }
                            }
                            if (!closeCoreChecked) {
                                if (!d.closeCoreMaps_halfR.isEmpty()) {
                                    result.put(d.closeCoreMaps_halfR.get(0), new ArrayList<>());
                                    result.get(d.closeCoreMaps_halfR.get(0)).add(d);
                                } else {
                                    result.put(d.closeCoreMaps_R.get(0), new ArrayList<>());
                                    result.get(d.closeCoreMaps_R.get(0)).add(d);
                                }
                            }
                        }
                    }
                }
            }
        }
        return result;
    }
    
    private void scanForCore(Set<CorePoint> cores) {
        for (CorePoint c : cores) {
            for (Integer sIndex : all_slides.keySet()) {
                if (!c.scannedSlide.contains(sIndex)) {
                    ArrayList<C_Data> neighbor_half_R = c.closeNeighbors_halfR.get(sIndex);
                    if (neighbor_half_R == null) {
                        neighbor_half_R = new ArrayList<>();
                        c.closeNeighbors_halfR.put(sIndex, neighbor_half_R);
                    }
                    ArrayList<C_Data> neighbor_R = c.closeNeighbors_R.get(sIndex);
                    if (neighbor_R == null) {
                        neighbor_R = new ArrayList<>();
                        c.closeNeighbors_R.put(sIndex, neighbor_R);
                    }
                    ArrayList<C_Data> neighbor_3half_R = c.closeNeighbors_3halfR.get(sIndex);
                    if (neighbor_3half_R == null) {
                        neighbor_3half_R = new ArrayList<>();
                        c.closeNeighbors_3halfR.put(sIndex, neighbor_3half_R);
                    }
                    ArrayList<C_Data> neighbor_2R = c.closeNeighbors_2R.get(sIndex);
                    if (neighbor_2R == null) {
                        neighbor_2R = new ArrayList<>();
                        c.closeNeighbors_2R.put(sIndex, neighbor_2R);
                    }
                    
                    for (C_Data d : all_slides.get(sIndex)) {
                        if (!d.closeCoreMaps_halfR.contains(c)
                                && !d.closeCoreMaps_R.contains(c)) {
                            double distance = DistanceFunction.euclideanDistance(c, d);
                            
                            if (distance <= Constants.R / 2) {
                                d.closeCoreMaps_halfR.add(c);
                                c.neighborCount += 1;
                                
                                c.closeNeighbors_halfR.get(sIndex).add(d);
                                
                            } else if (distance <= Constants.R) {
                                d.closeCoreMaps_R.add(c);
                                c.neighborCount += 1;
                                
                                c.closeNeighbors_R.get(sIndex).add(d);

//                break;
                            } else if (distance <= Constants.R * 1.5) {
//                

                                c.closeNeighbors_3halfR.get(sIndex).add(d);
                                
                            } else if (distance <= Constants.R * 2) {
                                
                                c.closeNeighbors_2R.get(sIndex).add(d);
                                
                            }
                            
                        }
                    }
                    
                    c.scannedSlide.add(sIndex);
                }
            }
        }
    }
    
    private void probe_slide(C_Data d, CorePoint c,  int slideIndex) {
        
        int countNeighbor = 0;

        //scan possible points 
        ArrayList<ArrayList<C_Data>> possibleCandidates = new ArrayList<>();
        if (!d.closeCoreMaps_halfR.isEmpty()) { // c is close neigbor
            assert d.closeCoreMaps_halfR.get(0)==c;
            
            //grab close neighbor in range R/2 of c
            
            countNeighbor += c.closeNeighbors_halfR.get(slideIndex).size();
            
            if (slideIndex < d.sIndex) {
                
                d.pred_neighbor_count.put(slideIndex, countNeighbor);
                if (neighborCountTrigger.containsKey(slideIndex)) {
                    neighborCountTrigger.get(slideIndex).add(d);
                } else {
                    HashSet<C_Data> hs = new HashSet<>();
                    hs.add(d);
                    neighborCountTrigger.put(slideIndex, hs);
                }
            } else {
                d.numSucceedingNeighbor += countNeighbor;
            }
            
            if (d.numSucceedingNeighbor + countNeighbor >= Constants.k) {
                return;
            }
            
            possibleCandidates.add(c.closeNeighbors_R.get(slideIndex));
            possibleCandidates.add(c.closeNeighbors_3halfR.get(slideIndex));
            if (c.closeNeighbors_R.get(slideIndex) == null) {
                System.out.println("WTF!!!!!");
            }
            if (c.closeNeighbors_3halfR.get(slideIndex) == null) {
                System.out.println("WTF!!!!!");
            }
        } else {
            assert c==d.closeCoreMaps_R.get(0);
            
            possibleCandidates.add(c.closeNeighbors_halfR.get(slideIndex));
            possibleCandidates.add(c.closeNeighbors_R.get(slideIndex));
            possibleCandidates.add(c.closeNeighbors_3halfR.get(slideIndex));
            possibleCandidates.add(c.closeNeighbors_2R.get(slideIndex));
        }
        
        long start = Utils.getCPUTime();
        for (ArrayList<C_Data> ps : possibleCandidates) {
//            System.out.println("Slide Index = "+ slideIndex);
//            System.out.println("Possibl candidate sizes = "+ ps.size());
            if (ps == null) {
                System.out.println("Slide index = " + slideIndex);
                System.out.println("Possible candidate size = " + possibleCandidates.size());
                
            }
            for (C_Data d2 : ps) {
                if (check_distance_neighbor_boolean(d, d2)) {
                    countNeighbor += 1;
                    if (d2.sIndex >= d.sIndex) {
                        d.numSucceedingNeighbor += 1;
                    }
                    if (countNeighbor >= Constants.k || d.numSucceedingNeighbor >= Constants.k) {
                        
                        break;
                    }
                }
            }
            if (countNeighbor >= Constants.k) {
                
                break;
            }
        }
        timeCheckingCandidates += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        d.neighborCount += countNeighbor;
        if (slideIndex < d.sIndex) {
//            d.numSucceedingNeighbor += countNeighbor;

            d.pred_neighbor_count.put(slideIndex, countNeighbor);
//            start = Utils.getCPUTime();
            if (neighborCountTrigger.containsKey(slideIndex)) {
                neighborCountTrigger.get(slideIndex).add(d);
            } else {
                HashSet<C_Data> hs = new HashSet<>();
                hs.add(d);
                neighborCountTrigger.put(slideIndex, hs);
            }
//            timeAddingToNeighborCount += (Utils.getCPUTime() - start) * 1.0 / 1000000000;
        }
//        System.out.println("Time looping through candidates = "+ ((Utils.getCPUTime() - start) * 1.0 / 1000000000));

    }
    
    private void probe(C_Data d, CorePoint c, int newestSlide) {
        if (d.lastProbRight < newestSlide) {
            //prob right first
            int slideIndex = d.lastProbRight + 1;
            if (d.lastProbRight == -1) {
                slideIndex = d.sIndex;
                
            }
            
            while (slideIndex <= newestSlide && d.neighborCount < Constants.k) {
                probe_slide(d, c, slideIndex);
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
                probe_slide(d, c, slideIndex);
                d.lastProbLeft = slideIndex;
                slideIndex--;
            }
        }
        
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
    }

    private HashMap<CorePoint, ArrayList<C_Data>> notifyCloseCore() {
        HashMap<CorePoint, ArrayList<C_Data>> result = new HashMap<>();

        for (Integer sIndex : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIndex)) {
                if (!(d instanceof CorePoint)) {
                    if (d.neighborCount < Constants.k) {
                        boolean needCheck = true;
                        if (!d.closeCoreMaps_halfR.isEmpty()) {
                            for (CorePoint c : d.closeCoreMaps_halfR) {
                                if (c.totalCloseNeighbor() >= Constants.k) {
                                    needCheck = false;
                                }
                            }
                        }
                        
                        if (needCheck) {
                            boolean closeCoreChecked = false;
                            for (CorePoint c : d.closeCoreMaps_halfR) {
                                if (result.containsKey(c)) {
                                    closeCoreChecked = true;
                                    result.get(c).add(d);
                                    break;
                                }
                                
                            }
                            if (!closeCoreChecked) {
                                for (CorePoint c : d.closeCoreMaps_R) {
                                    if (result.containsKey(c)) {
                                        closeCoreChecked = true;
                                        result.get(c).add(d);
                                        break;
                                    }
                                }
                            }
                            if (!closeCoreChecked) {
                                if (!d.closeCoreMaps_halfR.isEmpty()) {
                                    result.put(d.closeCoreMaps_halfR.get(0), new ArrayList<>());
                                    result.get(d.closeCoreMaps_halfR.get(0)).add(d);
                                } else {
                                    result.put(d.closeCoreMaps_R.get(0), new ArrayList<>());
                                    result.get(d.closeCoreMaps_R.get(0)).add(d);
                                }
                            }
                        }
                    }
                }
            }
        }
        return result;
    }

//    private HashMap<CorePoint,ArrayList<C_Data>>  findPointNeedProbe(HashMap<CorePoint, ArrayList<C_Data>> core_c_data_map) {
//        HashMap<CorePoint,ArrayList<C_Data>> results = new HashMap<>();
//        
//        for(CorePoint c: core_c_data_map.keySet()){
//            if(c.totalCloseNeighbor() < Constants.k){
//                
//            }
//        }
//        for (Integer sIndex : all_slides.keySet()) {
//            for (C_Data d : all_slides.get(sIndex)) {
//                //check if d is close neighbor of a full core
//                if (d.closeCoreMaps_halfR.isEmpty() || (d.closeCoreMaps_halfR.get(0).closeNeighbors_halfR.size() < Constants.k)) {
//                    if (d.neighborCount < Constants.k) {
//                        results.add(d);
//                    }
//                }
//            }
//        }
//        return results;
//    }
    final class ResultFindCore {
        
        private final double distance;
        private final ArrayList<CorePoint> cores;
        
        public ResultFindCore(double distance, ArrayList<CorePoint> cores) {
            this.distance = distance;
            this.cores = cores;
        }
        
        public double getDistance() {
            return this.distance;
        }
        
        public ArrayList<CorePoint> getCore() {
            return this.cores;
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
        
        private ArrayList<CorePoint> closeCoreMaps_halfR = new ArrayList<>();
        private ArrayList<CorePoint> closeCoreMaps_R = new ArrayList<>();
        private ArrayList<CorePoint> closeCoreMaps_3halfR = new ArrayList<>();
        private ArrayList<CorePoint> closeCoreMaps_2R = new ArrayList<>();
//        private CorePoint closeCore;
//        public ArrayList<CorePoint> linkCoreMaps = new ArrayList<>();

        public ArrayList<Integer[]> groups = new ArrayList<>();
        
        public double[] mean_list;
        public double[] std_list;
        public ArrayList<Integer> ordered_attributes = new ArrayList<>();
        
        public int sIndex = -1;
        
        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;
            
            this.sIndex = (arrivalTime - 1) / Constants.slide;
            
            if (null != Constants.dataFile) {
                switch (Constants.dataFile) {
                    
                    case "household2.txt":
                        
                        if (Constants.useLB1 || Constants.useUB1) {
                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6});
                        }
                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(2, 3, 6, 0, 5, 1, 4));
                        }
                        break;
                    case "new_hpc.txt":
                        this.groups.add(new Integer[]{0, 1, 3});
                        this.groups.add(new Integer[]{2, 4, 5, 6});
                        break;
                    case "covtype.data":
                        if (Constants.useLB1 || Constants.useUB1) {
//                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//                                21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
//                                41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54});
//                            this.groups.add(new Integer[]{ });
                            this.groups.add(new Integer[]{0});
                            this.groups.add(new Integer[]{1});
                            this.groups.add(new Integer[]{2});
                            this.groups.add(new Integer[]{4});
                            this.groups.add(new Integer[]{3});
                            this.groups.add(new Integer[]{5});
                            this.groups.add(new Integer[]{6});
                            this.groups.add(new Integer[]{7});
                            this.groups.add(new Integer[]{8});
                            this.groups.add(new Integer[]{9});
                            this.groups.add(new Integer[]{11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                                41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54});
//                            this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5});
//                            this.groups.add(new Integer[]{6, 7, 8, 9, 54});
                        }
                        
                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                                    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54));
                        }
                        break;
                    case "ethylene.txt":
                        
                        if (Constants.useLB1 || Constants.useUB1) {
                            
                            this.groups.add(new Integer[]{0, 2});
//                    this.groups.add(new Integer[]{1});
                            this.groups.add(new Integer[]{3, 4, 6, 7});
                            this.groups.add(new Integer[]{8, 14, 15});
                            this.groups.add(new Integer[]{1, 9});
                            this.groups.add(new Integer[]{10, 11, 12});
                            this.groups.add(new Integer[]{13});
                            this.groups.add(new Integer[]{5});
                        }
                        
                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7));
                        }
                        break;
                    
                    case "new_tao.txt":
                        this.groups.add(new Integer[]{0, 1, 2});
                        
                        break;
                    
                    case "tao.txt":
                        if (Constants.useLB1 || Constants.useUB1) {
                            this.groups.add(new Integer[]{0, 1, 2});
                        }
                        if (Constants.useLB2) {
                            this.ordered_attributes = new ArrayList<>(Arrays.asList(1, 2));
                        }
                        break;
                    
                    case "gaussian.txt":
                        if (Constants.useLB1 || Constants.useUB1) {
                            this.groups.add(new Integer[0]);
                        }
                    default:
                        break;
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
                        if (this.groups.get(i).length == 1) {
                            this.std_list[i] = 0;
                        } else {
                            temp = 0;
                            for (Integer j : this.groups.get(i)) {
                                temp += (this.values[j] - this.mean_list[i]) * (this.values[j] - this.mean_list[i]);
                            }
                            
                            this.std_list[i] = Math.sqrt(temp / this.groups.get(i).length);
                        }
                    }
                }
            }

//            this.isOutlier = true;
        }
        
        public C_Data() {
            
        }
        
        public int countNeighbor() {
            return neighborCount;
        }
        
    }
    
    class CorePoint extends C_Data {
        
        public HashSet<Integer> scannedSlide = new HashSet<>();
        
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_halfR = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_R = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_3halfR = new HashMap<>();
        public HashMap<Integer, ArrayList<C_Data>> closeNeighbors_2R = new HashMap<>();
        
        public int getTotalCoverPoint() {
            int t = 0;
            for (Integer sIndex : all_slides.keySet()) {
                if (closeNeighbors_halfR.containsKey(sIndex)) {
                    t += closeNeighbors_halfR.get(sIndex).size();
                }
                if (closeNeighbors_R.containsKey(sIndex)) {
                    t += closeNeighbors_R.get(sIndex).size();
                }
                if (closeNeighbors_3halfR.containsKey(sIndex)) {
                    t += closeNeighbors_3halfR.get(sIndex).size();
                }
                if (closeNeighbors_2R.containsKey(sIndex)) {
                    t += closeNeighbors_2R.get(sIndex).size();
                }
                
            }
            return t;
        }
//  
//        private int numPointsCovered = 0;

        public CorePoint(C_Data d) {
            this.sIndex = d.sIndex;
            this.values = d.values;
            this.hashCode = d.hashCode;
            this.arrivalTime = d.arrivalTime;
            this.groups = d.groups;
            this.mean_list = d.mean_list;
            this.std_list = d.std_list;
        }
        
        private int totalCloseNeighbor() {
            int t = 0;
            for (Integer sIndex : this.closeNeighbors_halfR.keySet()) {
                t += this.closeNeighbors_halfR.get(sIndex).size();
            }
            return t;
        }
        
    }
    
}
