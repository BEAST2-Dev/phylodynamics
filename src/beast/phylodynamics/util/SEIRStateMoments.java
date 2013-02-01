package beast.phylodynamics.util;

/**
 * Methods for calculating basic moments from lists of SIR trajectories
 *
 * @author tvaughan
 *
 */

import beast.phylodynamics.epidemiology.SEIRState;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

public class SEIRStateMoments {


    public static List<List<SEIRState>> getSample(List<List<SEIRState>> trajectoryList, double maxTime,  int Nsamples) throws Exception{

        List<List<SEIRState>> sampleList = new ArrayList<List<SEIRState>>();
        List<SEIRState> sample;
        double dt = maxTime / Nsamples;

        for (List<SEIRState> traj : trajectoryList){

            sample = new ArrayList<SEIRState>();
            ListIterator trajIt = traj.listIterator();
            SEIRState former = ((SEIRState) trajIt.next()).copy();// trajIt.remove();
            SEIRState current = ((SEIRState) trajIt.next()).copy(); //trajIt.remove();

            sample.add(former);

            int index = 0;
            while ( current.time < dt && trajIt.hasNext() ){
                index = trajIt.previousIndex();
                current = (SEIRState) trajIt.next();  //trajIt.remove();
            }

            while (traj.get(index).time > 0)
                if (traj.get(index).time > former.time)
                    index--;
                else break;

            former = traj.get(index).copy();

            for (int i = 1; i < Nsamples; i++) {

                if (former.time <= (dt*i) && current.time >= (dt*i)){
                    if (former.time == (dt*i))
                        sample.add(former);
                    else if (current.time == (dt*i))
                        sample.add(current);
                    else
                        sample.add(getmeanState(former, current, dt*i));
                }

                if (!trajIt.hasNext()) break;

                index = trajIt.previousIndex();
                current = (SEIRState) trajIt.next(); //trajIt.remove();
                
                while (current.time < dt*(i+1) && trajIt.hasNext() ){
                    index = trajIt.previousIndex();
                    current = ((SEIRState) trajIt.next()).copy();   //trajIt.remove();
                }
                while (traj.get(index).time > dt*i)
                    if (traj.get(index).time > former.time)
                        index--;
                    else break;

                
                former = traj.get(index).copy();
                
            }
            

            SEIRState last = traj.get(traj.size()-1);
            last.time = sample.get(sample.size()-1).time +dt;
            sample.add(last);

            sampleList.add(sample);

        }

        return sampleList;
    }

    protected static SEIRState getmeanState(SEIRState former, SEIRState current, double time){

        double c = (time-former.time) / (current.time - former.time);

        return new SEIRState((int)Math.round((1-c)*former.S+c*current.S),(int)Math.round((1-c)*former.E+c*current.E), (int)Math.round((1-c)*former.I+c*current.I), (int)Math.round((1-c)*former.R+c*current.R), time);

    }


    /**
     * Calculate mean of trajectories
     *
     * @param trajectoryList	list of trajectories returned by HybridTauleapSIR::genTrajectory
     * @param Nsamples          max number of samples in trajectory list
     * @return list containing the mean value at each time step
     */
    public static List<SEIRStateDouble> getMeans(List<List<SEIRState>> trajectoryList, int Nsamples) {

        List<SEIRStateDouble> means = new ArrayList<SEIRStateDouble>();

        int j;

        int Ntraj = trajectoryList.size();
//        int Nt = trajectoryList.get(0).size();

        Boolean timeSet;

        for (int i=0; i<Nsamples; i++) {

            SEIRStateDouble thismean = new SEIRStateDouble(0,0,0,0,0.0);
            timeSet = false;
            for (int traj=0; traj<Ntraj; traj++) { 

                j = i<trajectoryList.get(traj).size()?i:trajectoryList.get(traj).size()-1;

                thismean.S += trajectoryList.get(traj).get(j).S;
                thismean.E += trajectoryList.get(traj).get(j).E;
                thismean.I += trajectoryList.get(traj).get(j).I;
                thismean.R += trajectoryList.get(traj).get(j).R;

                if (!timeSet && i<trajectoryList.get(traj).size()) {
                    thismean.time = trajectoryList.get(traj).get(i).time;
                    timeSet = true;

                }
            }

            thismean.S /= Ntraj;
            thismean.E /= Ntraj;
            thismean.I /= Ntraj;
            thismean.R /= Ntraj;

            means.add(thismean);

        }

        return means;

    }

    /**
     * Calculate variance of trajectories
     *
     * @param trajectoryList	list of trajectories returned by HybridTauleapSIR::genTrajectory
     * @param means				list returned by SIRStateMoments::getMeans
     * @param Nsamples       max number of samples in trajectory list
     * @return list containing the variance at each time step
     */
    public static List<SEIRStateDouble> getVariances(List<List<SEIRState>> trajectoryList,
                                                     List<SEIRStateDouble> means, int Nsamples) {

        List<SEIRStateDouble> vars = new ArrayList<SEIRStateDouble>();

        int j;
        
        int Ntraj = trajectoryList.size();
//        int Nt = trajectoryList.get(0).size();

        Boolean[] timeSet = new Boolean[Nsamples];

        for (int i=0; i<Nsamples; i++) {

            SEIRStateDouble thisvar = new SEIRStateDouble(0,0,0,0,0.0);
            timeSet[i] = false;
            for (int traj=0; traj<Ntraj; traj++) {

                j = i<trajectoryList.get(traj).size()?i:trajectoryList.get(traj).size()-1;

                thisvar.S += trajectoryList.get(traj).get(j).S*trajectoryList.get(traj).get(j).S;
                thisvar.E += trajectoryList.get(traj).get(j).E*trajectoryList.get(traj).get(j).E;
                thisvar.I += trajectoryList.get(traj).get(j).I*trajectoryList.get(traj).get(j).I;
                thisvar.R += trajectoryList.get(traj).get(j).R*trajectoryList.get(traj).get(j).R;


                if (!timeSet[i] && i<trajectoryList.get(traj).size()) {
                    thisvar.time = trajectoryList.get(traj).get(i).time;
                    timeSet[i] = true;
                }

            }

            thisvar.S /= Ntraj;
            thisvar.E /= Ntraj;
            thisvar.I /= Ntraj;
            thisvar.R /= Ntraj;

            thisvar.S -= means.get(i).S*means.get(i).S;
            thisvar.E -= means.get(i).E*means.get(i).E;
            thisvar.I -= means.get(i).I*means.get(i).I;
            thisvar.R -= means.get(i).R*means.get(i).R;

//            thisvar.time = trajectoryList.get(0).get(i).time;


            vars.add(thisvar);

        }

        return vars;

    }

}
