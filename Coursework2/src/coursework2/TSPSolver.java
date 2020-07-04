package coursework2;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class TSPSolver {

    static boolean flag;

    public static ArrayList<City> readFile(String filename) {
        ArrayList<City> cities = new ArrayList<>();
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line = null;
            while ((line = in.readLine()) != null) {
                String[] blocks = line.trim().split("\\s+");
                if (blocks.length == 3) {
                    City c = new City();
                    c.city = Integer.parseInt(blocks[0]);
                    c.x = Double.parseDouble(blocks[1]);
                    c.y = Double.parseDouble(blocks[2]);
                    //System.out.printf("City %s %f %f\n", c.city, c.x, c.y);
                    cities.add(c);
                } else {
                    continue;
                }
            }
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
        }
        City.distances = new double[cities.size()][cities.size()];
        for (int i = 0; i < cities.size(); i++) {
            City ci = cities.get(i);
            for (int j = i; j < cities.size(); j++) {
                City cj = cities.get(j);
                City.distances[i][j] = City.distances[j][i] = Math.sqrt(Math.pow((ci.x - cj.x), 2) + Math.pow((ci.y - cj.y), 2));
            }
        }
        return cities;
    }

    public static ArrayList<City> solveProblem(ArrayList<City> citiesToVisit) {
        ArrayList<City> routine = new ArrayList<City>();
        City start = null;
        City current = null;
        // get city 0;
        for (int i = 0; i < citiesToVisit.size(); i++) {
            if (citiesToVisit.get(i).city == 0) {
                start = current = citiesToVisit.remove(i);
                routine.add(current);
                break;
            }
        }
        if (current == null) {
            System.out.println("Your problem instance is incorrect! Exiting...");
            System.exit(0);
        }
        // visit cities
        while (!citiesToVisit.isEmpty()) {
            double minDist = Double.MAX_VALUE;
            int index = -1;
            for (int i = 0; i < citiesToVisit.size(); i++) {
                double distI = current.distance(citiesToVisit.get(i));
                // index == -1 is needed in case the distance is really Double.MAX_VALUE.
                if (index == -1 || distI < minDist) {
                    index = i;
                    minDist = distI;
                }
            }
            //int index = 0;

            current = citiesToVisit.remove(index);
            routine.add(current);
        }
        routine.add(start); // go back to 0
        return routine;
    }

    public static double printSolution(ArrayList<City> routine) {
        double totalDistance = 0.0;
        for (int i = 0; i < routine.size(); i++) {
            if (i != routine.size() - 1) {
                System.out.print(routine.get(i).city + "->");
                totalDistance += routine.get(i).distance(routine.get(i + 1));
            } else {
                System.out.println(routine.get(i).city);
            }
        }
        return totalDistance;
    }

    /*
        Just evaluate the total distance. A simplified version of printSolution()
     */
    public static double evaluateRoutine(ArrayList<City> routine) {
        double totalDistance = 0.0;
        for (int i = 0; i < routine.size() - 1; i++) {
            totalDistance += routine.get(i).distance(routine.get(i + 1));
        }
        return totalDistance;
    }

    /*
        Moves the city at index "from" to index "to" inside the routine
     */
    private static void moveCity(ArrayList<City> routine, int from, int to) {
        // provide your code here.
        if (to - from != 1 && from != to) {
            City c = new City();
            c = routine.get(from);
            routine.remove(from);
            if (from < to) {
                routine.add(to - 1, c);
            } else {
                routine.add(to, c);
            }
        }
    }

    public static double evalMove(ArrayList<City> routine, int from, int to) {
        double olddisstance = routine.get(from - 1).distance(routine.get(from)) + routine.get(from).distance(routine.get(from + 1)) + routine.get(to - 1).distance(routine.get(to));
        double newdisstance = routine.get(from - 1).distance(routine.get(from + 1)) + routine.get(to - 1).distance(routine.get(from)) + routine.get(from).distance(routine.get(to));
        return olddisstance - newdisstance;

    }

    public static boolean moveFirstImprove(ArrayList<City> routine) {
        for (int i = 1; i < routine.size() - 1; i++) {
            for (int j = 1; j < routine.size(); j++) {
                if (i == j || i + 1 == j) {
                    continue;
                }
                double diff = evalMove(routine, i, j);
                if (diff - 0.00001 > 0) { // I really mean diff > 0 here
                    moveCity(routine, i, j);
                    return true;
                }
            }
        }
        return false;
    }

    public static void swapCity(ArrayList<City> routine, int index1, int index2) {
        City c = new City();
        c = routine.get(index1);
        routine.set(index1, routine.get(index2));
        routine.set(index2, c);
    }

    public static double evalSwap(ArrayList<City> routine, int index1, int index2) {
        City c = new City();
        int beforeC1, afterC1, beforeC2, afterC2, C1, C2;
        C1 = routine.get(index1).city;
        C2 = routine.get(index2).city;
        beforeC1 = routine.get(index1 - 1).city;
        afterC1 = routine.get(index1 + 1).city;
        beforeC2 = routine.get(index2 - 1).city;
        afterC2 = routine.get(index2 + 1).city;
        if (index1 + 1 == index2) {
            double oldDistance = c.distances[beforeC1][C1] + c.distances[C2][afterC2];
            double newDistance = c.distances[beforeC1][C2] + c.distances[C1][afterC2];
            return oldDistance - newDistance;
        } else {
            double oldDistance = c.distances[beforeC1][C1] + c.distances[C1][afterC1] + c.distances[beforeC2][C2] + c.distances[C2][afterC2];
            double newDistance = c.distances[beforeC1][C2] + c.distances[C2][afterC1] + c.distances[beforeC2][C1] + c.distances[C1][afterC2];
            return oldDistance - newDistance;
        }
    }

    public static boolean swapFirstImprove(ArrayList<City> routine) {
        for (int i = 1; i < routine.size() - 1; i++) {
            for (int j = i + 1; j < routine.size() - 1; j++) {
                double diff = evalSwap(routine, i, j);
                if (diff - 0.00001 > 0) { // I really mean diff > 0 here
                    swapCity(routine, i, j);
                    return true;
                }
            }
        }
        return false;
    }

    //  reverse the city between index "start+1" to index "end-1" inside the routine
    public static void reverseCity(ArrayList<City> routine, int start, int end) {
        ArrayList<City> temp = new ArrayList<>();
        int s = Math.min(start, end);
        int e = Math.max(start, end);
        if (start != end) {
            for (int i = s + 1; i <= e - 1; i++) {
                temp.add(routine.get(i));
            }
            int j = e - 1;
            for (int i = 0; i < e - s - 1; i++) {
                routine.set(j, temp.get(i));
                j--;
            }
        }
    }

    // evaluate the distances after reverse
    public static double evalReverse(ArrayList<City> routine, int s, int e) {
        City c = new City();
        int beforeC1, afterC2, C1, C2;
        C1 = routine.get(s + 1).city;
        C2 = routine.get(e - 1).city;
        beforeC1 = routine.get(s).city;
        afterC2 = routine.get(e).city;

        double oldDistance = c.distances[beforeC1][C1] + c.distances[C2][afterC2];
        double newDistance = c.distances[beforeC1][C2] + c.distances[C1][afterC2];
        return oldDistance - newDistance;

    }

    // traversed all the cities until the reversal was found.
    public static boolean reverseFirstImprove(ArrayList<City> routine) {
        for (int i = 0; i <= routine.size() - 1; i++) {
            for (int j = i + 1; j <= routine.size() - 1; j++) {
                double diff = evalReverse(routine, i, j);
                if (diff - 0.00001 > 0 && j - i != 1 && j - i != 2) { // I really mean diff > 0 here
                    reverseCity(routine, i, j);
                    return true;
                }
            }
        }
        return false;
    }

    // Further optimize TSP based on greedy algorithm
    public static ArrayList<City> improveRoutine(ArrayList<City> routine) {
        long startTime0 = System.currentTimeMillis();
        long endTime = startTime0;
        int MaxTime = 290;
        // Max_Segment refers to the maximum number of segments in the shuffle() can be divided 
        int Max_Segment = 10;
        //Calculate the coefficients of the quadratic function
        double k = (1 - Max_Segment) / (double) (MaxTime * MaxTime);
        while (true) {
            ArrayList<City> temp = new ArrayList<City>(routine);
            temp = shuffle(temp, (double) (endTime - startTime0) / 1000, MaxTime, Max_Segment, k);
            //If the path is changed by shuffling
            if (flag == true) {
                while (reverseFirstImprove(temp));
            }
            //If the current solution is less than the global optimal solution, save it as the global optimal
            if (evaluateRoutine(temp) < evaluateRoutine(routine)) {
                routine = (ArrayList<City>) temp.clone();
            }
            endTime = System.currentTimeMillis();

            //If the running time exceeds 4 minutes and 50 seconds, break and return
            if ((endTime - startTime0) > MaxTime * 1000) {
                break;
            }
        }
        return routine;
    }

    public static ArrayList<City> shuffle(ArrayList<City> routine, double CurrTime, int MaxTime, int MaxSegment, double k) {
        flag = true;
        //The size of the fluctuation allowed to be optimized
        double rate = 0.06;
        double t = evaluateRoutine(routine) * rate;
        
        int len = routine.size();
        //Calculate the number of segments
        int segment = MaxSegment + (int) (k * CurrTime * CurrTime);
        int[] position = new int[segment + 1];
        Random r = new Random();

        position[0] = 1;
        //Randomly divide the routine into many segments
        for (int i = 1; i < position.length - 1; i++) {

            position[i] = (int) (r.nextInt(len / (segment - 1) - 2) + len / (segment - 1) * (i - 1) + 1);

        }
        position[segment] = len - 1;
        ArrayList<Integer> order = new ArrayList<Integer>();
        for (int i = 0; i < segment; i++) {
            order.add(i);
        }
        ArrayList<City> shuffledRoutine = new ArrayList<>();
        //Reorder these segments
        shuffledRoutine.add(routine.get(0));
        Collections.shuffle(order);
        //Assign the reordered routine to shuffledRoutine
        for (int i = 0; i < position.length - 1; i++) {
            ArrayList<City> temp = new ArrayList<>(routine.subList(position[order.get(i)], position[order.get(i) + 1]));
            shuffledRoutine.addAll(temp);
        }
        shuffledRoutine.add(routine.get(0));
        //If the routine cannot pass the filter after the change
        if (Math.exp((evaluateRoutine(routine) - evaluateRoutine(shuffledRoutine)) / t) < r.nextDouble()) {
            flag = false;
            return routine;
        }
        return shuffledRoutine;
    }
}
