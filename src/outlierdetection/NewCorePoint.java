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
import java.util.Random;
import mtree.tests.Data;
import mtree.utils.Constants;
import mtree.utils.Utils;

/**
 *
 * @author luan
 */
public class NewCorePoint {

    public static HashMap<Integer, ArrayList<C_Data>> all_slides = new HashMap<>();

    public static ArrayList<CorePoint> all_core_points = new ArrayList<>();

    public static HashMap<Integer, ArrayList<C_Data>> trigger_list = new HashMap<>();

    public static int currentTime;

    public static int expiredSlideIndex = -1;

    public static double max_percentage_core = 0.01;

    public ArrayList<Data> detectOutlier(ArrayList<Data> data, int _currentTime, int W, int slide) {

        ArrayList<Data> result = new ArrayList<>((int) (data.size() * 1.0 / 100));
        currentTime = _currentTime;
        expiredSlideIndex = (currentTime - 1) / slide - Constants.W / Constants.slide;
        System.out.println("Expired slide index = " + expiredSlideIndex);
        processExpiredData();

        for (int i = 0; i < data.size(); i++) {
            Data o = data.get(i);
            C_Data d = new C_Data(o);

            ArrayList<C_Data> idx = all_slides.get(d.sIndex);

            if (idx != null) {
                idx.add(d);
            } else {
                idx = new ArrayList<>(Constants.slide);
                idx.add(d);
                all_slides.put(d.sIndex, idx);
            }
        }

        int newestSlide = (currentTime - 1) / Constants.slide;
        System.out.println("Newest slide = " + newestSlide);

        if (data.size() == Constants.W) {
//            selectCore2();
//
//            scanForCore(all_core_points);
            scanAndFormCore();
            System.out.println("Finished selecting and scanning for core");
            System.out.println("Num core = " + all_core_points.size());
//            System.out.println("Filtering cores....");
//            all_core_points = filterCore(100);

//            for (Integer sIndex : all_slides.keySet()) {
//                ArrayList<CorePoint> newCores = promoteNewCore(sIndex);
//                scanForCore(newCores);
//
//            }
//            System.out.println("Num core point = " + all_core_points.size());
        } else if (data.size() == Constants.slide) {

            long start = Utils.getCPUTime();
            matchCore(newestSlide);
            System.out.println("Time for matching core = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
//            HashSet<CorePoint> coreNeedScan = findCoreToScan(newestSlide);
//            System.out.println("Number of core need to scan = " + coreNeedScan.size());
//            updateCore(coreNeedScan, newestSlide);
//            start = Utils.getCPUTime();
//            ArrayList<CorePoint> newCores = promoteNewCore(newestSlide);
//            System.out.println("Num New Cores = " + newCores.size());
//            System.out.println("Time for promoting new cores = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
//
//            start = Utils.getCPUTime();
//            scanForCore(newCores);
//
//            System.out.println("Time for scanning for new cores = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
//
//            start = Utils.getCPUTime();
//            matchCore(newestSlide, newCores);
//            System.out.println("Time for matching new cores = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

            start = Utils.getCPUTime();
            removeUselessCore();
            System.out.println("Time for removing useless cores = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

            start = Utils.getCPUTime();
            System.out.println("Num cores = " + all_core_points.size());
            populateStatus();
            System.out.println("Finished updating core");

            System.out.println("Time for populating status = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);

        }

//        System.out.println("Core Point Check");
//        for (CorePoint c : all_core_points) {
//            System.out.println("Total cover #points = " + c.getTotalCoverPoints());
//        }
//        System.out.println("------------------");
        long start = Utils.getCPUTime();
        int count = 0;
        for (Integer sIdx : all_slides.keySet()) {
//            System.out.println("Checking sIdx = " + sIdx);
            for (C_Data d : all_slides.get(sIdx)) {
                if (d.isOutlier) {
                    probNeighbors(d, newestSlide);
                    count += 1;
                    if (d.countNeighbor() < Constants.k) {
                        result.add(d);
                    }

                }

            }
        }
        System.out.println("Time for probing = " + (Utils.getCPUTime() - start) * 1.0 / 1000000000);
        System.out.println("# Probing Points = " + count);

        return result;
    }

//    private boolean isValidCore(C_Data d: )
    private void scanAndFormCore() {
        for (Integer sIndex : all_slides.keySet()) {
            for (C_Data d : all_slides.get(sIndex)) {
                if (d.closeCore == null) {
                    //promote d to be a core
                    CorePoint c = new CorePoint(d.values);
                    all_core_points.add(c);
                    d.closeCore = c;
                    HashSet<C_Data> closeNeighbor = new HashSet<>();
                    closeNeighbor.add(d);
                    c.closeNeighbors.put(d.sIndex, closeNeighbor);

                    for (Integer s2 : all_slides.keySet()) {
                        for (C_Data d2 : all_slides.get(s2)) {
                            if (d2.closeCore == null) {
                                double distance = DistanceFunction.euclideanDistance(c, d2);
                                if (distance <= Constants.R / 2) {
                                    d2.closeCore = c;
                                    closeNeighbor = c.closeNeighbors.get(d2.sIndex);
                                    if (closeNeighbor == null) {
                                        closeNeighbor = new HashSet<>();
                                        closeNeighbor.add(d2);
                                        c.closeNeighbors.put(d2.sIndex, closeNeighbor);
                                    } else {
                                        closeNeighbor.add(d2);
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }
    }

    private void removeUselessCore() {
        int count = 0;
        for (int i = all_core_points.size() - 1; i >= 0; i--) {
            if (all_core_points.get(i).getTotalCoverPoints() == 0) {
                all_core_points.remove(i);
                count += 1;
            }
        }
        System.out.println("Remove #cores = " + count);
    }

    private void populateStatus() {
        for (CorePoint c : all_core_points) {
//            int numCloseNeighbors = c.numPointsCovered;
            boolean status = true;
            if (c.numPointsCovered > Constants.k) {
                status = false;
            }
            for (Integer sIndex : c.closeNeighbors.keySet()) {
                for (C_Data d : c.closeNeighbors.get(sIndex)) {
                    d.isOutlier = status;

                }
            }
        }
    }

    private void processExpiredData() {

        all_slides.remove(expiredSlideIndex);
        for (CorePoint c : all_core_points) {
            c.closeNeighbors.remove(expiredSlideIndex);
//            c.doubleRNeighbors.remove(expiredSlideIndex);
        }

        if (trigger_list.containsKey(expiredSlideIndex)) {
            for (C_Data d : trigger_list.get(expiredSlideIndex)) {
                d.neighborCount -= d.pred_neighbor_count.get(expiredSlideIndex);
                d.pred_neighbor_count.remove(expiredSlideIndex);

            }
        }
//        for (Integer sIdx : all_slides.keySet()) {
//            for (C_Data d : all_slides.get(sIdx)) {
//                if (d.pred_neighbor_count.containsKey(expiredSlideIndex)) {
//                    d.pred_neighbor_count.remove(expiredSlideIndex);
//                }
//            }
//        }

    }

    public static List<C_Data> randomSubList(ArrayList<C_Data> list, int newSize) {

        Collections.shuffle(list);
        return list.subList(0, newSize);
    }

    public static double[] computeAverageValue(List<C_Data> selectedData) {
        int dim = selectedData.get(0).dimensions();
        double[] result = new double[dim];
        for (C_Data d : selectedData) {
            for (int i = 0; i < dim; i++) {
                result[i] += d.values[i];
            }
        }
        for (int i = 0; i < dim; i++) {
            result[i] /= selectedData.size();
        }
        return result;
    }

    private ArrayList<CorePoint> promoteNewCore(int newestSlide) {
        ArrayList<CorePoint> results = new ArrayList<>();

        ArrayList<C_Data> candidates = new ArrayList<>();
        for (C_Data d : all_slides.get(newestSlide)) {
            if (d.closeCore == null) {
                candidates.add(d);
            }
        }
        for (C_Data d : candidates) {
            boolean isGood = true;
            for (CorePoint c : all_core_points) {
                double distance = DistanceFunction.euclideanDistance(d, c);
                if (distance <= Constants.R) {
                    isGood = false;
                }
            }
            if (isGood) {
                CorePoint c = new CorePoint(d.values);
                all_core_points.add(c);
                results.add(c);
                if (all_core_points.size() > 500) {
                    break;
                }
            }
        }

        return results;
    }

    private void selectCore2() {
        ArrayList<Integer> all_slide_incides = new ArrayList<>(all_slides.keySet());
        Collections.sort(all_slide_incides, Collections.reverseOrder());
        int max_num_core = 50;
        Random r = new Random();
        while (all_core_points.size() < max_num_core) {
            int slide_idx = r.nextInt(all_slide_incides.size());
            int selectedSlide = all_slide_incides.get(slide_idx);
            int selectedOffset = r.nextInt(Constants.slide);
            C_Data selectedPoint = all_slides.get(selectedSlide).get(selectedOffset);
            boolean isGoodCore = true;
            for (CorePoint c : all_core_points) {
                Double distance = DistanceFunction.euclideanDistance(selectedPoint, c);
                if (distance <= Constants.R) {
                    isGoodCore = false;
                    break;
                }
            }
            if (isGoodCore) {
                all_core_points.add(new CorePoint(selectedPoint.values));
                System.out.println("Found " + all_core_points.size());
            }

        }
    }

    private void selectCore() {
        ArrayList<Integer> all_slide_incides = new ArrayList<>(all_slides.keySet());
        Collections.sort(all_slide_incides, Collections.reverseOrder());

        int max_num_core = (int) (Constants.W * max_percentage_core);

        int cur_num_core = 0;

        Random r = new Random();
        while (cur_num_core < max_num_core) {
            int slide_idx = r.nextInt(all_slide_incides.size());
            List<C_Data> selectedData = randomSubList(all_slides.get(slide_idx),
                    1);
            double[] averageValue = computeAverageValue(selectedData);
            CorePoint c = new CorePoint(averageValue);
            all_core_points.add(c);
            cur_num_core++;

        }

    }

    public void probNeighbors(C_Data d, int newestSlide) {

        //populate neighbors
//        d.neighborCount = 0;
//        for (int i = newestSlide - Constants.W / Constants.slide + 1; i <= newestSlide; i++) {
//            if (d.lastProbLeft == -1 || i < d.lastProbLeft) {
//                if (d.closeCore != null && d.closeCore.closeNeighbors.containsKey(i)) {
//                    if (d.closeCore.closeNeighbors.containsKey(i)) {
//                        d.neighborCount += d.closeCore.closeNeighbors.get(i).size();
//                    }
//                }
//            } else if (i >= d.lastProbLeft && i < d.sIndex) {
//                d.neighborCount += d.pred_neighbor_count.get(i);
//            } else if (i >= d.sIndex && i <= d.lastProbRight) {
//                d.neighborCount += d.numSucceedingNeighbor;
//            } else if (i > d.lastProbRight) {
//                if (d.closeCore!=null && d.closeCore.closeNeighbors.containsKey(i)) {
//                    d.neighborCount += d.closeCore.closeNeighbors.get(i).size();
//                }
//            }
//        }
        //prob right first 
        int idx = d.lastProbRight + 1;
        if (idx == 0) {
            idx = d.sIndex;//first time probing
        }
        while (idx <= newestSlide && d.countNeighbor() < Constants.k) {
            int countNeighbor = prob_slide(d, idx);

            d.lastProbRight = idx;
            d.numSucceedingNeighbor += countNeighbor;
            if (d.closeCore != null && d.closeCore.closeNeighbors.containsKey(idx)) {
                d.neighborCount += countNeighbor - d.closeCore.closeNeighbors.get(idx).size();
            } else {
                d.neighborCount += countNeighbor;
            }
            if (countNeighbor >= Constants.k) {
                return;
            }
            idx++;
        }
        //prob left
        idx = d.lastProbLeft - 1;
        if (idx < -1) {
            idx = d.sIndex - 1;
        }
        while (idx > expiredSlideIndex && d.countNeighbor() < Constants.k) {
            int countNeighbor = prob_slide(d, idx);
            d.lastProbLeft = idx;
            d.pred_neighbor_count.put(idx, countNeighbor);
            if (d.closeCore != null && d.closeCore.closeNeighbors.containsKey(idx)) {
                d.neighborCount += countNeighbor - d.closeCore.closeNeighbors.get(idx).size();
            } else {
                d.neighborCount += countNeighbor;
            }

            //put d to trigger list
            if (trigger_list.containsKey(idx)) {
                trigger_list.get(idx).add(d);
            } else {
                ArrayList<C_Data> points = new ArrayList<>();
                points.add(d);
                trigger_list.put(idx, points);
            }

            if (countNeighbor >= Constants.k) {
                return;
            }
            idx--;
        }

    }

    public int prob_slide(C_Data d, int sIndex) {
        int countNeighbor = 0;
        //stop when finding enough k neighbors
        CorePoint selectedCore = null;
//        int bestCloseNeighbor = -1;

        //select the core with the largest close neighbors in the slide
//            for (CorePoint c : d.closeCoreMaps) {
//                if (c.closeNeighbors.containsKey(sIndex)
//                        && !c.closeNeighbors.get(sIndex).isEmpty()) {
//                    if (c.closeNeighbors.get(sIndex).size() > bestCloseNeighbor) {
//                        bestCloseNeighbor = c.closeNeighbors.get(sIndex).size();
//                        selectedCore = c;
//                    }
//                    if (c.closeNeighbors.get(sIndex).size() >= Constants.k) {
////                    System.out.println("Found #close neighbors =  " + c.closeNeighbors.get(sIndex).size());
//                        return c.closeNeighbors.get(sIndex).size();
//                    }
//                }
//            }
//            System.out.println("Found #close neighbors =  " + selectedCore.closeNeighbors.get(sIndex).size());
//            if (selectedCore != null) {
        //probe points using double R candidates of selected Core 
        if (d.closeCore.closeNeighbors.containsKey(sIndex)) {
            selectedCore = d.closeCore;
            HashSet<C_Data> probed = new HashSet<>(selectedCore.closeNeighbors.get(sIndex));
            countNeighbor = selectedCore.closeNeighbors.get(sIndex).size();
            ArrayList<C_Data> candidates = all_slides.get(sIndex);
            for (C_Data d2 : candidates) {
                if (!probed.contains(d2)) {
                    double distance = DistanceFunction.euclideanDistance(d, d2);
                    if (distance <= Constants.R) {
                        countNeighbor += 1;
                        if (countNeighbor >= Constants.k) {
                            return countNeighbor;
                        }
                    }

                }
            }
        } else {
            ArrayList<C_Data> candidates = all_slides.get(sIndex);
            for (C_Data d2 : candidates) {

                double distance = DistanceFunction.euclideanDistance(d, d2);
                if (distance <= Constants.R) {
                    countNeighbor += 1;
                    if (countNeighbor >= Constants.k) {
                        return countNeighbor;
                    }
                }

            }
        }
//            }

        //        if (selectedCore == null && !d.linkCoreMaps.isEmpty() && d.linkCoreMaps.get(0).checkedSlide.contains(sIndex)) {
        //            //select the first core to probe doubleRcandidates
        //            selectedCore = d.linkCoreMaps.get(0);
        //            ArrayList<C_Data> candidates = selectedCore.doubleRNeighbors.get(sIndex);
        //            for (C_Data d2 : candidates) {
        //                double distance = DistanceFunction.euclideanDistance(d, d2);
        //                if (distance <= Constants.R) {
        //                    countNeighbor += 1;
        //                    if (countNeighbor >= Constants.k) {
        //                        return countNeighbor;
        //                    }
        //                }
        //
        //            }
        //        } 
//        else if (selectedCore == null) {
//            ArrayList<C_Data> candidates = all_slides.get(sIndex);
//            for (C_Data d2 : candidates) {
//                double distance = DistanceFunction.euclideanDistance(d, d2);
//                if (distance <= Constants.R) {
//                    countNeighbor += 1;
//                    if (countNeighbor >= Constants.k) {
//                        return countNeighbor;
//                    }
//                }
//
//            }
//        }
        return countNeighbor;
    }

    private void scanForCore(ArrayList<CorePoint> corepoints) {
        for (CorePoint c : corepoints) {
            for (Integer sIdx : all_slides.keySet()) {
//                ArrayList<C_Data> doubleRNeighbors = new ArrayList<>();
                HashSet<C_Data> closeNeighbors = new HashSet<>();
                for (C_Data d : all_slides.get(sIdx)) {
//                    if (d.closeCore == null) {
                    double distance = DistanceFunction.euclideanDistance(c, d);
//                    if (distance <= Constants.R * 2) {
//                        doubleRNeighbors.add(d);
//
//                        if (distance <= Constants.R) {
//                            d.linkCoreMaps.add(c);
                    if (distance <= Constants.R / 2) {
                        closeNeighbors.add(d);

                        d.closeCore = c;
                    }
//                    }
//                        }
//
//                    }
                }
                c.closeNeighbors.put(sIdx, closeNeighbors);
//                c.doubleRNeighbors.put(sIdx, doubleRNeighbors);
                c.checkedSlide.add(sIdx);
            }
        }
    }

    private void updateCore(HashSet<CorePoint> coreList, int newestSlide) {
        for (CorePoint c : coreList) {

//            ArrayList<C_Data> doubleRNeighbors = new ArrayList<>();
            HashSet<C_Data> closeNeighbors = c.closeNeighbors.get(newestSlide);
            if (closeNeighbors == null) {
                closeNeighbors = new HashSet<>();
            }

            for (C_Data d : all_slides.get(newestSlide)) {
                if (!closeNeighbors.contains(d)) {
                    double distance = DistanceFunction.euclideanDistance(c, d);
//                    if (distance <= Constants.R * 2) {
//                        doubleRNeighbors.add(d);
                    if (distance <= Constants.R / 2) {
                        closeNeighbors.add(d);

                        d.closeCore = c;

                        if (closeNeighbors.size() > Constants.k) {
                            break;
                        }
                    }
//                        else if (distance <= Constants.R) {
//                            d.linkCoreMaps.add(c);
//
//                        }

//                    }
                }
            }
            c.closeNeighbors.put(newestSlide, closeNeighbors);
//            c.doubleRNeighbors.put(newestSlide, doubleRNeighbors);
            c.checkedSlide.add(newestSlide);

        }
    }

    private ArrayList<CorePoint> filterCore(int remainCore) {
        ArrayList<ArrayList<C_Data>> pointsCoveredByCore = new ArrayList<>();
        for (CorePoint c : all_core_points) {
            pointsCoveredByCore.add(c.getCoveredPoints());
        }
        ArrayList<CorePoint> selectedCores = new ArrayList<>();
//        int selected = 0;
        HashSet<C_Data> coveredPoints = new HashSet<>();
//        HashSet<CorePoint> checked = new HashSet<>();
        while (true) {
            int bestNewlyCovered = 0;
            CorePoint bestCore = null;
            for (int i = 0; i < all_core_points.size(); i++) {
                CorePoint c = all_core_points.get(i);
                if (!selectedCores.contains(c)) {
                    HashSet<C_Data> clone_coveredPoints = new HashSet<>();
                    clone_coveredPoints = (HashSet) coveredPoints.clone();
                    //try to add c
                    clone_coveredPoints.addAll(pointsCoveredByCore.get(i));
                    int newCovered = clone_coveredPoints.size() - coveredPoints.size();
                    if (newCovered > bestNewlyCovered) {
                        bestNewlyCovered = newCovered;
                        bestCore = c;
                    }
                }
            }
            if (bestNewlyCovered > 0) {
                coveredPoints.addAll(bestCore.getCoveredPoints());
                selectedCores.add(bestCore);
                System.out.println("Selected " + selectedCores.size());
                if (selectedCores.size() >= remainCore) {
                    break;
                }
            } else {
                break;
            }
        }
        return selectedCores;
    }

    private void matchCore(int newestSlide, ArrayList<CorePoint> cores) {
        for (C_Data d : all_slides.get(newestSlide)) {
            if (d.closeCore == null) {
                for (CorePoint c : cores) {
                    double distance = DistanceFunction.euclideanDistance(d, c);
                    if (distance <= Constants.R / 2) {
                        d.closeCore = c;
                        HashSet<C_Data> closeNeighbors = c.closeNeighbors.get(newestSlide);
                        if (closeNeighbors == null) {
                            closeNeighbors = new HashSet<>();
                            closeNeighbors.add(d);
                            c.closeNeighbors.put(newestSlide, closeNeighbors);
                        } else {
                            closeNeighbors.add(d);
                        }
                        break;
                    }
                }

            }
        }
    }

    private void matchCore(int newestSlide) {

        for (C_Data d : all_slides.get(newestSlide)) {
            if (d.closeCore == null) //find close core for each d
            {
                for (CorePoint c : all_core_points) {

                    double distance = DistanceFunction.euclideanDistance(d, c);
                    if (distance <= Constants.R / 2) {
                        d.closeCore = c;
                        HashSet<C_Data> closeNeighbors = c.closeNeighbors.get(newestSlide);
                        if (closeNeighbors == null) {
                            closeNeighbors = new HashSet<>();
                            closeNeighbors.add(d);
                            c.closeNeighbors.put(newestSlide, closeNeighbors);
                        } else {
                            closeNeighbors.add(d);
                        }
                        break;
                    }
//                else if (distance <= Constants.R) {
//                    d.linkCoreMaps.add(c);
//                }
                }
                //if cannot find any core
                if (d.closeCore == null) {
                    CorePoint c = new CorePoint(d.values);
                    all_core_points.add(c);
                    d.closeCore = c;
                    HashSet<C_Data> closeNeighbors = new HashSet<>();
                    closeNeighbors.add(d);
                    c.closeNeighbors.put(newestSlide, closeNeighbors);

                    closeNeighbors.add(d);

                }
            }
        }

    }

    private HashSet<CorePoint> findCoreToScan(int newestSlide) {

        HashSet<CorePoint> results = new HashSet<>();
        for (C_Data d : all_slides.get(newestSlide)) {
            if (d.closeCore != null) {
                if (d.closeCore.getTotalCoverPoints() <= Constants.k
                        && d.closeCore.closeNeighbors.get(newestSlide).size() < Constants.k) {
                    results.add(d.closeCore);
                }
            }
        }
//        for (CorePoint c : all_core_points) {
//
//            if (c.closeNeighbors.containsKey(newestSlide)
//                    && !c.closeNeighbors.get(newestSlide).isEmpty()
//                    && c.closeNeighbors.get(newestSlide).size() <= Constants.k) {
//                results.add(c);
//            }
//        }
        return results;
    }

    class C_Data extends Data {

        private int numSucceedingNeighbor = 0;
        private boolean isOutlier;
        private int lastProbRight = -1;
        private int lastProbLeft = -1;
        private HashMap<Integer, Integer> pred_neighbor_count = new HashMap<>();
        private int neighborCount = 0;

//        private ArrayList<CorePoint> closeCoreMaps = new ArrayList<>();
        private CorePoint closeCore;
//        public ArrayList<CorePoint> linkCoreMaps = new ArrayList<>();

        public int sIndex = -1;

        public C_Data(Data d) {
            this.arrivalTime = d.arrivalTime;
            this.values = d.values;
            this.hashCode = d.hashCode;

            this.sIndex = (arrivalTime - 1) / Constants.slide;
            this.isOutlier = true;
        }

        public C_Data() {

        }

        public int countNeighbor() {
//            int totalNeighbor = numSucceedingNeighbor;
//            for (Integer sIdx : pred_neighbor_count.keySet()) {
//                totalNeighbor += pred_neighbor_count.get(sIdx);
//            }
//            return totalNeighbor;
            return neighborCount;
        }

    }

    class CorePoint extends C_Data {

        public HashMap<Integer, HashSet<C_Data>> closeNeighbors = new HashMap<>();
//        public HashMap<Integer, ArrayList<C_Data>> doubleRNeighbors = new HashMap<>();

        public ArrayList<Integer> checkedSlide = new ArrayList<>();

        private int numPointsCovered = 0;

        public CorePoint(double... values) {
            this.values = values;
        }

        public ArrayList<C_Data> getCoveredPoints() {
            ArrayList<C_Data> result = new ArrayList<>();
            for (Integer sIndex : closeNeighbors.keySet()) {
                result.addAll(closeNeighbors.get(sIndex));
            }
            return result;
        }

        public int getTotalCoverPoints() {
            int result = 0;
            for (Integer sIndex : closeNeighbors.keySet()) {
                result += closeNeighbors.get(sIndex).size();
            }
            numPointsCovered = result;
            return result;
        }
    }

}
