/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outlierdetection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import mtree.tests.Data;
import mtree.tests.MesureMemoryThread;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class UpperLower {

    class ULObject extends Data {

        public boolean isOutlier;
        public boolean isSafe;
        public int numSucEvidence;
        public int numPreEvidence = -1;
        public double mean;
        public double std;
        public double[] a;
        public int last_probed_slide = -1;
        public int last_max_probed_slide = -1;
        public int numbins = 10;
        public ArrayList[] distance_bins = new ArrayList[numbins];
        
        

//        public double[]
        public ULObject(Data d, int currentTime) {
            super();
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
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
            this.a = new double[this.values.length];
            for (int i = 0; i < this.values.length; i++) {
                this.a[i] = values[i] - mean;
            }

            numPreEvidence = -1;

            for (int i = 0; i < distance_bins.length; i++) {
                distance_bins[i] = new ArrayList<>();
            }
        }

//        public int getCurrentSlideIndex() {
//            return (int) Math.floor((arrivalTime - 1) / Constants.slide);
//        }
        public void updateSuccEvidence() {
            numSucEvidence++;
        }

        public boolean isInlier() {

            return numSucEvidence + numPreEvidence >= Constants.k;

        }

    }
    
    public static double MSU_j(double mx, double my, double sx, double sy, double[] a, double[] b, int j, int d){
        return 0;
    }
    
    public static double MSL_j(double mx, double my, double sx, double sy, double[] a, double[] b, int j, int d){
        double distance = 0;
        double epsilon = 0.000000001;
        distance += d * ((mx - my) * (mx - my) + (sx - sy) * (sx - sy));
        
        for (int i = 0; i < j; i++){
            distance += sx * sy * (b[i] / (sy + epsilon) - a[i] / (sx + epsilon)) * (b[i] / (sy + epsilon) - a[i] / (sx + epsilon));
        }
        return Math.sqrt(distance);
        
    }

    public static double MSU(double mx, double my, double sx, double sy, double[] a, double[] b, int j, int d, double prev_distance) {
        double distance = 0;
        double epsilon = 0.000000001;
        if (j == 0) {
            distance += d * ((mx - my) * (mx - my) + (sx + sy) * (sx + sy));
        } else {
            distance = prev_distance - sx * sy * (b[j] / (sy + epsilon) + a[j] / (sx + epsilon)) * (b[j] / (sy + epsilon) + a[j] / (sx + epsilon));
        }
        return Math.sqrt(distance);
    }

    public static double MSL(double mx, double my, double sx, double sy, double[] a, double[] b, int j, int d, double prev_distance) {
        double distance = 0;
        if (j == 0) {
            distance += d * ((mx - my) * (mx - my) + (sx - sy) * (sx - sy));
        } else {
            distance = prev_distance + sx * sy * (b[j] / sy - a[j] / sx) * (b[j] / sy - a[j] / sx);
        }
        return Math.sqrt(distance);
    }

    public static double MSU(ULObject x, ULObject y, int j, double prev_distance) {
//        System.out.println("x.mean = " + x.mean);
//        System.out.println("y.mean = " + y.mean);
//        System.out.println("x.std = " + x.std);
//        System.out.println("y.std = " + y.std);
//        System.out.println("Dimensions = "+ x.dimensions());
        return MSU(x.mean, y.mean, x.std, y.std, x.a, y.a, j, x.dimensions(), prev_distance);
    }

    public static double MSL(ULObject x, ULObject y, int j, double prev_distance) {
        return MSL(x.mean, y.mean, x.std, y.std, x.a, y.a, j, x.dimensions(), prev_distance);
    }

    class Slide {

        public List<ULObject> points = new ArrayList<>();
        public int id;

        public Window window;
        public HashSet<ULObject> triggered = new HashSet<>();

        public Slide(ArrayList<Data> data, int currentTime) {

            int start_time = currentTime - data.size();
            for (int i = 0; i < data.size(); i++) {
                ULObject obj = new ULObject(data.get(i), start_time + i);
                points.add(obj);
            }

        }

    }

    class Window {

        public ArrayList<Slide> slides = new ArrayList<>();
        public int startSlide;

        public Slide getNewestSlide() {
            if (slides.isEmpty()) {
                return null;
            }
            return slides.get(slides.size() - 1);
        }

        public Window() {

            startSlide = 0;

        }

        public Slide getExpiringSlide() {
            if (slides.size() <= Constants.W / Constants.slide) {
                return null;
            } else {
                return slides.get(0);
            }
        }

//        public Slide getExpiredSlide() {
//
//            if (slides.size() <= Constants.W / Constants.slide) {
//                return null;
//            } else {
//                return slides.get(0);
//            }
//        }
        public void addNewSlide(Slide s) {
            Slide newestSlide = getNewestSlide();
            if (newestSlide != null) {
                s.id = newestSlide.id + 1;
            } else {
                s.id = 0;
            }
            s.window = this;
            slides.add(s);
            if (slides.size() >= Constants.W / Constants.slide) {
                startSlide++;
            }
        }

        public ArrayList<ULObject> findNeighbors(ULObject q) {
            q.isSafe = false;
            q.isOutlier = false;
            q.numPreEvidence = 0;
            q.numSucEvidence = 0;
//            HashSet<Integer> filter_list = new HashSet<>(Constants.W);
            int[] filter_list = new int[Constants.W];

            double[] all_upper_distances = new double[Constants.W];
            double[] all_lower_distances = new double[Constants.W];
            for (int i = 0; i < filter_list.length; i++) {
                filter_list[i] = 0;
                all_lower_distances[i] = 0;
                all_upper_distances[i] = 0;
            }
            int start_time = slides.get(0).points.get(0).arrivalTime;
//            HashMap<Integer, Double> all_upper_distances = new HashMap<>(Constants.W);
//            HashMap<Integer, Double> all_lower_distances = new HashMap<>(Constants.W);
//            int succeeding_neighbors = 0;
//            int preceding_neighbors = 0;
            q.numPreEvidence = 0;
            q.numSucEvidence = 0;
            int q_slide = q.arrivalTime / Constants.slide;
            ArrayList<ULObject> neighbors = new ArrayList<>();
            boolean canStop = false;
            for (int j = 0; j < q.dimensions(); j++) {
                for (int i = slides.size() - 1; i >= 0; i--) {
                    Slide s = slides.get(i);
                    for (int t = 0; t < s.points.size(); t++) {
                        ULObject obj = s.points.get(t);
                        if (obj.arrivalTime == q.arrivalTime) {
                            continue;
                        }
                        if (filter_list[obj.arrivalTime - start_time] == 1) {
                            continue;
                        }
//                        if (filter_list.contains(obj.arrivalTime - start_time)){
//                            continue;
//                        }
                        double msu_d = 0;
                        double msl_d = 0;
                        if (j == 0) {
                            msu_d = MSU(q, obj, 0, 0);
                            msl_d = MSL(q, obj, 0, 0);

                            all_upper_distances[obj.arrivalTime - start_time] = msu_d;
                            all_lower_distances[obj.arrivalTime - start_time] = msl_d;
                        } else {
//                            System.out.println("previous distance = " +  all_upper_distances.get(obj.arrivalTime));
                            msu_d = MSU(q, obj, j, all_upper_distances[obj.arrivalTime - start_time]);
                            msl_d = MSU(q, obj, j, all_lower_distances[obj.arrivalTime - start_time]);
                            all_upper_distances[obj.arrivalTime - start_time] = msu_d;
                            all_lower_distances[obj.arrivalTime - start_time] = msl_d;
                        }
                        //filtering using msu and msl
                        if (msu_d < Constants.R) {
                            neighbors.add(obj);
                            if (neighbors.size() >= Constants.k) {
                                canStop = true;

                            }

                            if (q_slide <= i + startSlide) {
                                q.numSucEvidence += 1;
                            } else {
                                q.numPreEvidence += 1;
                            }

                            //update last probe
                            if (i + this.startSlide < q.last_probed_slide || q.last_probed_slide == -1) {
                                q.last_probed_slide = i + this.startSlide;
                            }
                            filter_list[obj.arrivalTime - start_time] = 1;
//                            filter_list.add(obj.arrivalTime - start_time);

                            //put neighbor to distance_bins
                            int bin_idx = (int) (msu_d / Constants.R * q.numbins);

                            q.distance_bins[bin_idx].add(obj);
                            int count = 0;
                            //find neighbors from neighbors of obj
                            for (int bin_i = 0; bin_i < obj.numbins - bin_idx; bin_i++) {
                                ArrayList<ULObject> candidates = obj.distance_bins[bin_i];
                                System.out.println("Candidate size "+candidates.size());
                                for (ULObject u : candidates) {
                                    if (u.arrivalTime >= start_time && u.arrivalTime != q.arrivalTime && neighbors.contains(u) == false) {
                                        neighbors.add(u);
                                        q.distance_bins[bin_i+bin_idx].add(u);
                                        count +=1;
                                        if (u.arrivalTime/Constants.slide >= q.arrivalTime/Constants.slide) {
                                            q.numSucEvidence += 1;
                                        } else {
                                            q.numPreEvidence += 1;
                                        }
                                        if (neighbors.size() >= Constants.k) {
                                            canStop = true;

                                            //update last probe
                                            if (u.arrivalTime / Constants.slide < q.last_probed_slide || q.last_probed_slide == -1) {
                                                q.last_probed_slide = u.arrivalTime / Constants.slide;
                                            }
                                        }
                                    }
                                }
                            }
                            System.out.println("Found "+count +" neighbors!");
                        }
                        if (msl_d > Constants.R) {
                            filter_list[obj.arrivalTime - start_time] = 1;
//                            filter_list.add(obj.arrivalTime - start_time);
                        }

//                        if (canStop == true) {
//
//                            break;
//                        }
                    }
                    if (canStop == true) {

                        break;
                    }
                }
                if (canStop == true) {

                    break;
                }
            }

            q.isOutlier = neighbors.size() < Constants.k;
            q.isSafe = q.numSucEvidence >= Constants.k;
            if (q.isSafe == false && q.last_probed_slide >= this.startSlide) {
                this.slides.get(q.last_probed_slide - this.startSlide).triggered.add(q);
            }
            return neighbors;
        }

    }

    public ArrayList<Data> outlierList;
    public Window window = new Window();

    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int currentTime, int W, int slide) {

        //clear points in expiring slide
//        Slide expiringSlide = window.getExpiringSlide();
//        if (expiringSlide != null) {
//            expiringSlide.points.clear();
//            //Remove expiring slide from window
//            window.slides.remove(0);
//        }
        outlierList = new ArrayList<>();
        long startCPUTime = Utils.getCPUTime();

        if (data.size() == Constants.W) {
            // split into slides
            int numSlide = (int) Math.ceil(Constants.W * 1.0 / Constants.slide);
            for (int i = 0; i < numSlide; i++) {

                ArrayList<Data> d = new ArrayList<>();
                for (int j = 0; j < Constants.slide; j++) {
                    if (i * Constants.slide + j < data.size()) {
                        d.add(data.get(i * Constants.slide + j));
                    }
                }
                Slide s = new Slide(d, currentTime);
                window.addNewSlide(s);

            }

        } else if (data.size() <= Constants.slide) {
            // add this slide to window
            Slide s = new Slide(data, currentTime);
            window.addNewSlide(s);

        }
        long currentCPUTime = Utils.getCPUTime();
        MesureMemoryThread.timeForNewSlide += currentCPUTime - startCPUTime;

        //Remove expire slide
        Slide expiredSlide = window.getExpiringSlide();
        if (expiredSlide != null) {
            expiredSlide.points.clear();
            window.slides.remove(0);
            expiredSlide.triggered.forEach((obj) -> {

                //Reprob neighbors for obj
                if (!obj.isSafe) {
                    window.findNeighbors(obj);
                }
            });
        }

        if (currentTime == Constants.W) {
            //find neighbors of all data
            window.slides.forEach((s) -> {
                s.points.forEach((obj) -> {
                    window.findNeighbors(obj);
                });
            });
        } else if (currentTime > Constants.W) {
            //find neighbors for the new slide
            window.getNewestSlide().points.forEach((obj) -> {
                window.findNeighbors(obj);
            });
        }

        //run outlier detection
        window.slides.forEach((s) -> {
            s.points.stream().filter((obj) -> (obj.isOutlier)).forEachOrdered((obj) -> {
                outlierList.add(obj);
            });
        });

        return outlierList;

    }

}
