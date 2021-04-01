/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author Luan Tran
 */
public class MCOD_MESI_Upper_Lower {

//    public static HashMap<Cluster, HashSet<MCData>> associates = new HashMap<>();
    public static ArrayList<Cluster> microClusters = new ArrayList<>();
    public static PriorityQueue<MCData> event_queue = new PriorityQueue(new MCDataNeighborComparator());
    public static HashMap<Integer, ArrayList<MCData>> PDList = new HashMap<>();

    // public static int numberWindows = 0;
    public static int currentTime;

    public static ArrayList<MCData> dispersedList = new ArrayList<>();

    public static ArrayList<Double> p_add_to_cluster = new ArrayList<>();
    public static ArrayList<Double> p_distance_filtered_by_LB = new ArrayList<>();
    public static ArrayList<Double> p_distance_filered_by_UB = new ArrayList<>();
    public static long all_distance_computations = 0;
    public static long filtered_b_LB = 0;

    public static long all_distance_after_LB = 0;
    public static long filtered_by_UB = 0;

    public static long time_adding_to_cluster = 0;

    public static long time_finding_neighbors = 0;
    public static long count_data = 0;
//    public static ArrayList<Integer> number_points_in_clusters = new ArrayList<>();

    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {
        ArrayList<Data> result = new ArrayList<>();
        currentTime = _currentTime;

        if (data.size() == slide) {
            int expiredSlideIndex = (currentTime - Constants.W - 1) / slide;

            processExpiredData(expiredSlideIndex);
        }
//        int added_to_cluster = 0;
//        int newestSlideIdx = (currentTime - 1) / slide;
        ArrayList<MCData> to_process = new ArrayList<>();
        for (Data o : data) {
            MCData d = new MCData(o);
            to_process.add(d);
//            assert newestSlideIdx == d.getSlideIndex();
            ArrayList<MCData> newestSlide = PDList.get(d.getSlideIndex());
            if (newestSlide == null) {
                newestSlide = new ArrayList<>();
                newestSlide.add(d);
                PDList.put(d.getSlideIndex(), newestSlide);
            } else {
                newestSlide.add(d);
            }

        }

        for (MCData d : to_process) {
            processData(d, false);

//            if (!d.clusters.isEmpty()) {
//                added_to_cluster += 1;
//            }
        }

//        if (data.size() > 0) {
//            p_add_to_cluster.add(added_to_cluster * 1.0 / data.size());
//        }
        //find outlier     
        for (Integer slideIndex : PDList.keySet()) {
            ArrayList<MCData> datas = PDList.get(slideIndex);
            int newestSlideIdx = (currentTime - 1) / Constants.slide;
//            System.out.println("Newest slide idx = " + newestSlideIdx);
            datas.stream().filter((d) -> (d.isOutlier())).map((d) -> {
                reProbe(d, newestSlideIdx);
                return d;
            }).filter((d) -> (d.isOutlier())).forEach((d) -> {
                result.add(d);
            });
        }
        dispersedList.clear();
        return result;

    }

    public void reProbe(MCData d, int newestSlideIndex) {
        int slideIndex = d.lastProbe + 1;
        while (slideIndex <= newestSlideIndex) {

            for (Cluster c : microClusters) {
                if ( //                        DistanceFunction.euclideanDistance(d, c.center) <= Constants.R * 3 / 2
                        check_distance_3R2_neighbor(d, c.center) <= Constants.R * 3 / 2) {
                    ArrayList<MCData> points = c.slide_members.get(slideIndex);
                    if (points != null) {
                        for (MCData d2 : points) {
                            if ( //                                    DistanceFunction.euclideanDistance(d, d2) <= Constants.R
                                    check_distance_neighbor(d, d2) <= Constants.R) {
                                d.numSucceedingNeighbors++;
                            }
                        }

                    }
                }
            }

            //probe in PD list
            ArrayList<MCData> points = PDList.get(slideIndex);
            for (MCData d2 : points) {
                if (check_distance_neighbor(d, d2) <= Constants.R) {
                    d.numSucceedingNeighbors++;
                }
            }

            if (!d.isOutlier()) {
                d.lastProbe = slideIndex;
                return;
            }
            d.lastProbe = slideIndex;
            slideIndex++;
        }
    }

    private Cluster findClusterToAdd(MCData d) {

        Cluster result = null;
//        double bestDistance = Double.MAX_VALUE;
        for (Cluster cluster : microClusters) {

//            double dis = DistanceFunction.euclideanDistance(d, cluster.center);
            double dis = check_distance_R2_neighbor(d, cluster.center);
            if (dis <= Constants.R / 2) {
//                    if (dis < bestDistance) {
//                        bestDistance = dis;
                result = cluster;
                break;
            }

        }
        return result;

    }

    private void processExpiredData(int expiredSlide) {
        if (PDList.get(expiredSlide) != null) {
            ArrayList<MCData> points = PDList.get(expiredSlide);
            for (MCData p : points) {
                p.clean();
                p.clusters.clear();
                p = null;
            }
            points.clear();
        }

        PDList.remove(expiredSlide);

        processEventQueue(expiredSlide);
        PDList.values().stream().forEach((datas) -> {
            datas.stream().forEach((d) -> {
                if (d.numPrecedingNeighbor.containsKey(expiredSlide)) {
                    HashSet<MCData> points = d.numPrecedingNeighbor.get(expiredSlide);
                    if (points != null) {
                        for (MCData p : points) {
                            p.clean();
                        }
                        points.clear();
                    }
                    d.numPrecedingNeighbor.remove(expiredSlide);
                }
            });
        });
        //update neighbor list
        dispersedList.clear();
        //remove from clusters 
        for (int i = microClusters.size() - 1; i >= 0; i--) {
            Cluster c = microClusters.get(i);
            //c.slide_members
            c.removeExpiredSlide(expiredSlide);
            //check if c has less than K+1 memebers 
            if (c.notEnoughMembers()) {
                //proces dispersed clusters
                for (ArrayList<MCData> members : c.slide_members.values()) {

                    members.stream().map((d) -> {
                        d.clusters.remove(c);
                        return d;
                    }).filter((d) -> (d.clusters.isEmpty())).map((d) -> {
                        d.clean();
                        return d;
                    }).forEach((d) -> {
                        dispersedList.add(d);
                    });

                }
                c.slide_members.clear();
                microClusters.remove(c);
            }
        }
        //process dispersed list
        dispersedList.stream().forEach((d) -> {
            processData(d, true);
        });
        dispersedList.clear();
        //remove from PD 

    }

    private void processEventQueue(int expiredSlide) {
        while (true) {
            if (event_queue.isEmpty()) {
                break;
            }
            MCData d = event_queue.peek();
            if (d.getOldestSlide() <= expiredSlide || d.getSlideIndex() <= expiredSlide
                    || d.numPrecedingNeighbor.isEmpty()) {
                event_queue.poll();
                if (d.numPrecedingNeighbor != null) {
                    //d.updateEarliestNeighbor();
                    for (int i = d.getOldestSlide(); i <= expiredSlide; i++) {
                        if (d.numPrecedingNeighbor.get(i) != null) {
                            d.numPrecedingNeighbor.get(i).clear();

                            d.numPrecedingNeighbor.remove(i);
                        }
                    }
//                    d.updateEarliestNeighbor();
                }
                if (d.isUnsafeInlier() && !d.numPrecedingNeighbor.isEmpty() && d.getSlideIndex() > expiredSlide) {
                    event_queue.add(d);
                } else if (d.getSlideIndex() <= expiredSlide) {
                    d.clean();
                }
            } else {
                break;
            }
        }
    }

    private void processData(MCData d, boolean isFromDispersedList) {
        long start = Utils.getCPUTime();
        //find cluster to add
        Cluster cluster = findClusterToAdd(d);
        time_adding_to_cluster += Utils.getCPUTime() - start;
        //add to cluster
        if (cluster != null) {
            start = Utils.getCPUTime();
            cluster.addNewPointToSlideMember(d);
            time_adding_to_cluster += Utils.getCPUTime() - start;
            PDList.get(d.getSlideIndex()).remove(d);

        } //add to PD list/event queue
        else {
            start = Utils.getCPUTime();
            //find neighbor for d 
            probe(d, isFromDispersedList);
            time_finding_neighbors += Utils.getCPUTime() - start;
            if (d.clusters.isEmpty()) {

                if (d.isUnsafeInlier()) {
                    event_queue.add(d);
                }
            }

        }
        count_data += 1;
    }

    public double check_distance_neighbor(MCData d, MCData d2) {

//        System.out.println("LB_D = " + DistanceFunction.LB_distance2(d, d2));
//        System.out.println("UB_D = " + DistanceFunction.UB_distance2(d, d2));
//        System.out.println("Exact d = " + DistanceFunction.euclideanDistance(d, d2));
//        System.out.println("*************************");
        if (Constants.useLB1) {
            double lb_d = DistanceFunction.LB_distance1(d, d2);
            if (lb_d > Constants.R ) {
                return lb_d;
            }
        } else if (Constants.useLB2) {
            double lb_d = DistanceFunction.LB_distance2(d, d2);
            if (lb_d > Constants.R) {
                return lb_d;
            }
        }

        if (Constants.useUB1) {
            double ub_d = DistanceFunction.UB_distance1(d, d2);
            if (ub_d <= Constants.R ) {
                return ub_d;
            }
        } else if (Constants.useUB2) {
            double ub_d = DistanceFunction.UB_distance2(d, d2);
            if (ub_d <= Constants.R ) {
                return ub_d;
            }
        }
//        
        double exact_d = DistanceFunction.euclideanDistance(d, d2);

        return exact_d;
    }

    public double check_distance_neighbor_reversed(MCData d, MCData d2) {

//        System.out.println("LB_D = " + DistanceFunction.LB_distance2(d, d2));
//        System.out.println("UB_D = " + DistanceFunction.UB_distance2(d, d2));
//        System.out.println("Exact d = " + DistanceFunction.euclideanDistance(d, d2));
//        System.out.println("*************************");
        double ub_d = DistanceFunction.UB_distance2(d, d2);
        if (ub_d <= Constants.R) {
            return ub_d;
        }

        double lb_d = DistanceFunction.LB_distance2(d, d2);
        if (lb_d > Constants.R) {
            return lb_d;
        }
//        
        double exact_d = DistanceFunction.euclideanDistance(d, d2);

        return exact_d;
    }

    public double check_distance_3R2_neighbor(MCData d, MCData d2) {
        
        if (Constants.useLB1) {
            double lb_d = DistanceFunction.LB_distance1(d, d2);
            if (lb_d > Constants.R * 3/2) {
                return lb_d;
            }
        } else if (Constants.useLB2) {
            double lb_d = DistanceFunction.LB_distance2(d, d2);
            if (lb_d > Constants.R * 3/2) {
                return lb_d;
            }
        }

        if (Constants.useUB1) {
            double ub_d = DistanceFunction.UB_distance1(d, d2);
            if (ub_d <= Constants.R * 3/2) {
                return ub_d;
            }
        } else if (Constants.useUB2) {
            double ub_d = DistanceFunction.UB_distance2(d, d2);
            if (ub_d <= Constants.R * 3/2) {
                return ub_d;
            }
        }
//        double lb_d = DistanceFunction.LB_distance2(d, d2);
//        if (lb_d > Constants.R * 3.0 / 2) {
//            return lb_d;
//        }
//        double ub_d = DistanceFunction.UB_distance2(d, d2);
//        if (ub_d <= Constants.R * 3.0 / 2) {
//            return ub_d;
//        }
        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public double check_distance_3R2_neighbor_reversed(MCData d, MCData d2) {

        double ub_d = DistanceFunction.UB_distance2(d, d2);
        if (ub_d <= Constants.R * 3.0 / 2) {
            return ub_d;
        }

        double lb_d = DistanceFunction.LB_distance2(d, d2);
        if (lb_d > Constants.R * 3.0 / 2) {
            return lb_d;
        }
        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public double check_distance_R2_neighbor(MCData d, MCData d2) {
        if (Constants.useLB1) {
            double lb_d = DistanceFunction.LB_distance1(d, d2);
            if (lb_d > Constants.R/ 2) {
                return lb_d;
            }
        } else if (Constants.useLB2) {
            double lb_d = DistanceFunction.LB_distance2(d, d2);
            if (lb_d > Constants.R /2) {
                return lb_d;
            }
        }

        if (Constants.useUB1) {
            double ub_d = DistanceFunction.UB_distance1(d, d2);
            if (ub_d <= Constants.R/2) {
                return ub_d;
            }
        } else if (Constants.useUB2) {
            double ub_d = DistanceFunction.UB_distance2(d, d2);
            if (ub_d <= Constants.R/2) {
                return ub_d;
            }
        }
        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public double check_distance_R2_neighbor_reversed(MCData d, MCData d2) {

        double ub_d = DistanceFunction.UB_distance2(d, d2);
        if (ub_d <= Constants.R * 1.0 / 2) {
            return ub_d;
        }

        double lb_d = DistanceFunction.LB_distance2(d, d2);
        if (lb_d > Constants.R * 1.0 / 2) {
            return lb_d;
        }
        double exact_d = DistanceFunction.euclideanDistance(d, d2);
        return exact_d;
    }

    public void probe(MCData d, boolean isFromDispersedList) {

        int slideIndex = d.getSlideIndex();
        ArrayList<MCData> closeNeighbors = new ArrayList<>();

        HashMap<Integer, HashSet<MCData>> slide_neighbors = new HashMap<>();

        //Form Cluster
        while ((slideIndex > d.getSlideIndex() - Constants.W / Constants.slide)
                && (closeNeighbors.size() < Constants.minSizeOfCluster - 1)
                && slideIndex >= 0) {
            HashSet<MCData> neighbors = new HashSet<>();
            //find neighbor in PD list
            ArrayList<MCData> datas = PDList.get(slideIndex);
//            System.out.println("Slide index = " + slideIndex);
//            System.out.println("d.slideIndex = " + d.getSlideIndex());
//            System.out.println("PDList key: ");
//            for (int k : PDList.keySet()) {
//                System.out.print(k + ", ");
//            }
//            System.out.println("");
            if (datas != null) {
                for (MCData d2 : datas) {
                    if (d2.arrivalTime != d.arrivalTime) {
                        double distance = check_distance_neighbor(d, d2);

                        if (distance <= Constants.R) {
                            neighbors.add(d2);

//                            distance = DistanceFunction.euclideanDistance(d, d2);
                            if (distance <= Constants.R / 2) {
                                closeNeighbors.add(d2);
                            }
                        }
                    }
                }
            }

            slide_neighbors.put(slideIndex, neighbors);

            if (closeNeighbors.size() >= Constants.minSizeOfCluster - 1) {
                //form new cluster
                formNewCluster(d, closeNeighbors);
                break;
            }
            slideIndex--;
        }

        if (d.clusters.isEmpty()) {

            slideIndex = d.getSlideIndex();
            int countNeighbor = 0;
            //find neighbors in clusters and join with neighbors in PD
            while ((slideIndex > d.getSlideIndex() - Constants.W / Constants.slide)
                    && (countNeighbor < Constants.k)) {

                HashSet<MCData> neighbors = slide_neighbors.get(slideIndex);
                //find neighbor in clusters 
                for (Cluster c : microClusters) {

                    if (check_distance_3R2_neighbor(d, c.center) <= Constants.R * 3 / 2 //                            DistanceFunction.euclideanDistance(d, c.center) <= Constants.R * 3 / 2
                            ) {
                        ArrayList<MCData> dataList = c.slide_members.get(slideIndex);
                        if (dataList != null) {
                            for (MCData d2 : dataList) {
                                if (d2.arrivalTime != d.arrivalTime) {
                                    double distance = check_distance_neighbor(d, d2);
                                    if (distance <= Constants.R) {
                                        neighbors.add(d2);
                                    }
                                }
                            }
                        }
                    }
                }

                if (neighbors != null) {
                    countNeighbor += neighbors.size();

                    if (slideIndex == d.getSlideIndex()) {
                        d.numSucceedingNeighbors += neighbors.size();
                    } else {
                        d.numPrecedingNeighbor.put(slideIndex, neighbors);
                    }

                    if (countNeighbor >= Constants.k) {

                        break;

                    }
                }
                slideIndex--;

            }
        }
    }

    private void formNewCluster(MCData d, ArrayList<MCData> closeNeighbors) {
        closeNeighbors.add(d);
        Cluster cluster = new Cluster();
        cluster.center = d;
        for (MCData d2 : closeNeighbors) {
            //remove from PD list
            PDList.get(d2.getSlideIndex()).remove(d2);
            if (d2.isUnsafeInlier() && d2.clusters.isEmpty()) {
                event_queue.remove(d2);
            }

            d2.clean();
            cluster.addNewPointToSlideMember(d2);
        }

        microClusters.add(cluster);
    }

//    private ArrayList<MCData> findNeighborInPD(MCData d, double R, int k) {
//
//        int slideIndex = d.getSlideIndex();
//        ArrayList<MCData> result = new ArrayList<>();
//        while (result.size() < k && slideIndex > d.getSlideIndex() - Constants.W / Constants.slide
//                && slideIndex >= 0) {
//            ArrayList<MCData> datas = PDList.get(slideIndex);
//            if (datas != null) {
//                datas.stream().filter((d2) -> (DistanceFunction.euclideanDistance(d, d2) <= R)).forEach((d2) -> {
//                    result.add(d2);
//                });
//            }
//            slideIndex--;
//        }
//        return result;
//    }
//    private ArrayList<MCData> findNeighborInCluster(MCData d, double R, int k, boolean useMaxCluster) {
//        ArrayList<MCData> result = new ArrayList<>();
//        int slideIndex = d.getSlideIndex();
//        while (slideIndex >= d.getSlideIndex() - Constants.W / Constants.slide
//                && slideIndex >= 0 && result.size() <= k) {
////            microClusters.stream().filter((c) -> (DistanceFunction.euclideanDistance(d, c.center) <= Constants.R / 2 + R)).map((c) -> c.slide_members.get(slideIndex)).filter((datas) -> (datas != null)).forEach((datas) -> {
////                datas.stream().filter((d2) -> (DistanceFunction.euclideanDistance(d, d2) <= R
////                        && (!useMaxCluster || d2.clusters.size() < Constants.maxClusterEachPoint))).forEach((d2) -> {
////                    result.add(d2);
////                });
////            });
//            for (Cluster c : microClusters) {
//                if (DistanceFunction.euclideanDistance(d, c.center) <= Constants.R / 2 + R) {
//                    for (MCData d2 : c.slide_members.get(slideIndex)) {
//                        if (DistanceFunction.euclideanDistance(d, d2) <= R
//                                && (!useMaxCluster || d2.clusters.size() < Constants.maxClusterEachPoint)) {
//                            result.add(d2);
//                        }
//                    }
//                }
//            }
//            slideIndex--;
//        }
//
//        return result;
//    }
    private static class MCDataNeighborComparator implements Comparator<MCData> {

        @Override
        public int compare(MCData t, MCData t1) {
            if (t.getOldestSlide() > t1.getOldestSlide()) {
                return 1;
            } else {
                return -1;
            }
        }
    }

    private class Cluster {

        public MCData center;
        public HashMap<Integer, ArrayList<MCData>> slide_members = new HashMap<>();

        public int numMembers() {
            int count = 0;
            count = slide_members.values().stream().map((members) -> members.size()).reduce(count, Integer::sum);
            return count;
        }

//        public boolean isSafe() {
//
//            return numMembers() > Constants.k;
//        }
        public Cluster() {
        }

        public void addNewPointToSlideMember(MCData d) {
            d.clusters.add(this);
            ArrayList<MCData> points = slide_members.get(d.getSlideIndex());
            if (points == null) {
                points = new ArrayList<>();
                points.add(d);
                slide_members.put(d.getSlideIndex(), points);
            } else {
                points.add(d);
            }

        }

        private void removeExpiredSlide(int expiredSlide) {
            if (this.slide_members.containsKey(expiredSlide)) {
                ArrayList<MCData> points = this.slide_members.get(expiredSlide);
                if (points != null) {
                    for (MCData p : points) {
                        p.clean();
                    }
                    points.clear();

                }

                this.slide_members.remove(expiredSlide);
            }

        }

        private int getTotalMembers() {
            int count = 0;
            count = slide_members.keySet().stream().map((key) -> slide_members.get(key).size()).reduce(count, Integer::sum);
            return count;
        }

        private boolean notEnoughMembers() {
            return this.getTotalMembers() < Constants.minSizeOfCluster;
        }
    }

    class MCData extends Data {

        public ArrayList<Cluster> clusters = new ArrayList<>();
        public int numSucceedingNeighbors;
        public HashMap<Integer, HashSet<MCData>> numPrecedingNeighbor;
        public int lastProbe;
        public int sIndex;

        public double[] a;
        public double mean;
        public double std;

        public ArrayList<Integer[]> groups = new ArrayList<>();

        public double[] mean_list;
        public double[] std_list;
        public int group_size = 2;
//        public int earliestExpireTime;

        public int getOldestSlide() {
            int result = getSlideIndex();
            for (Integer slide : numPrecedingNeighbor.keySet()) {
                if (slide < result) {
                    result = slide;
                }
            }
            return result;
        }

        public MCData(Data d) {
            super();
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;
            clusters = new ArrayList<>();
            numPrecedingNeighbor = new HashMap<>();
            numSucceedingNeighbors = 0;
            sIndex = (int) Math.floor((arrivalTime - 1) / Constants.slide);
            lastProbe = sIndex;
            this.mean = 0;
            for (int i = 0; i < this.values.length; i++) {
                this.mean += this.values[i];
            }
            this.mean = this.mean / this.values.length;

            this.std = 0;
            for (int i = 0; i < this.values.length; i++) {
                this.std += (this.values[i] - this.mean) * (this.values[i] - this.mean);
            }
            this.std = this.std / this.values.length;
            this.std = Math.sqrt(this.std);

//            this.a = new double[this.values.length];
//            for (int i = 0; i < this.values.length; i++) {
//                this.a[i] = this.values[i] - this.mean;
//            }
            if ("household2.txt".equals(Constants.dataFile)) {

                this.groups.add(new Integer[]{0, 5});
                this.groups.add(new Integer[]{1, 4});
                this.groups.add(new Integer[]{2});
                this.groups.add(new Integer[]{3, 6});
            } else if ("covtype.data".equals(Constants.dataFile)) {
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
                this.groups.add(new Integer[]{54});
                this.groups.add(new Integer[]{10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53});

            } else if ("ethylene.txt".equals(Constants.dataFile)) {
//                this.groups.add(new Integer[]{0, 2});
//                this.groups.add(new Integer[]{1});
//                this.groups.add(new Integer[]{3, 4, 6, 7});
//                this.groups.add(new Integer[]{8, 14, 15});
//                this.groups.add(new Integer[]{9});
//                this.groups.add(new Integer[]{10, 11, 12});
//                this.groups.add(new Integer[]{13});
                if (Constants.useLB1 || Constants.useUB1) {

//                        this.groups.add(new Integer[]{0, 2});
////                    this.groups.add(new Integer[]{1});
//                        this.groups.add(new Integer[]{3, 4, 6, 7});
//                        this.groups.add(new Integer[]{8, 14, 15});
//                        this.groups.add(new Integer[]{1, 9});
//                        this.groups.add(new Integer[]{10, 11, 12});
//                        this.groups.add(new Integer[]{13});
                    this.groups.add(new Integer[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
                        10, 11, 12, 13});
                }
//                break;

            } else if ("new_tao.txt".equals(Constants.dataFile)) {
                this.groups.add(new Integer[]{0, 1, 2});
            }
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

        public int countNeighborToSlide(int slideIndex) {
            int result = numSucceedingNeighbors;
            for (Integer slide : numPrecedingNeighbor.keySet()) {
                HashSet<MCData> datas = numPrecedingNeighbor.get(slide);
                result = result + datas.size();
            }

            return result;
        }

        private boolean isNeighbor(MCData d2) {
            return DistanceFunction.euclideanDistance(this, d2) <= Constants.R && this.arrivalTime != d2.arrivalTime;
        }

        private void addPrecedingNeighbor(MCData d2) {

            if (this.isOutlier()) {
                int slideIndex = d2.getSlideIndex();
                if (numPrecedingNeighbor.get(slideIndex) != null) {
                    HashSet<MCData> datas = numPrecedingNeighbor.get(slideIndex);
                    datas.add(d2);
                    // numPrecedingNeighbor.put(slideIndex, );
                } else {
                    HashSet<MCData> datas = new HashSet<>();
                    datas.add(d2);
                    numPrecedingNeighbor.put(slideIndex, datas);
                }
            }

        }

        private int getSlideIndex() {
            return sIndex;
        }

        private boolean isOutlier() {

            int numNeighbor = numSucceedingNeighbors;
            for (int slideIndex : numPrecedingNeighbor.keySet()) {
                numNeighbor += numPrecedingNeighbor.get(slideIndex).size();
            }
//            numNeighbor = numPrecedingNeighbor.keySet().stream().map((slideIndex) -> numPrecedingNeighbor.get(slideIndex)).map((datas) -> datas.size()).reduce(numNeighbor, Integer::sum);
            return numNeighbor < Constants.k;
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
        private boolean isUnsafeInlier() {

            return !isOutlier() && numSucceedingNeighbors < Constants.k;
        }

        private void clean() {

            numSucceedingNeighbors = 0;
            numPrecedingNeighbor.clear();
//            clusters.clear();
        }

        private boolean isInCluster() {
            return clusters.isEmpty();
        }
    }
}
