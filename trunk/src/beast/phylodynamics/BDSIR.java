package beast.phylodynamics;

import beast.core.Input;
import beast.core.Description;
import beast.core.Valuable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.speciation.BirthDeathSkylineModel;

import java.util.Arrays;

/**
 * @author Denise Kuhnert
 *
 */

@Description("Adaptation of Tanja Stadler's BirthDeathSamplingModel, to allow for birth and death rates to change at times t_i")
public class BDSIR extends BirthDeathSkylineModel {


    public Input<Valuable> S0_input =
            new Input<Valuable>("S0", "The numbers of susceptible individuals");

    public Input<RealParameter> m_dS =
            new Input<RealParameter>("dS", "dS vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);
    public Input<RealParameter> m_dE =
            new Input<RealParameter>("dE", "dE vector containing the changes in numbers of susceptibles per location");
    public Input<RealParameter> m_dR =
            new Input<RealParameter>("dR", "dR vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);

    public Input<Boolean> checkTreeConsistent = new Input<Boolean>("checkTreeConsistent", "check if trajectory is consistent with number of lineages in tree? default true", true);

    Double S0;
    Double[] dS;
    Double[] dE;
    Double[] dR;

    int dim;
    double T;
    int ntaxa;


    @Override
    public void initAndValidate() throws Exception {

        S0 =  (S0_input.get().getArrayValue());

        dS = m_dS.get().getValues();
        dim = dS.length;

        birthChanges = dim - 1;
        super.initAndValidate();

        if (transform){
            if (R0.get().getDimension() != 1)// || becomeUninfectiousRate.get().getDimension() != 1 || samplingProportion.get().getDimension() != 1)
                throw new RuntimeException("R0, becomeUninfectiousRate and samplingProportion have to be 1-dimensional!");
        } else {
            if (birthRate.get().getDimension() != 1)//  || death.length != 1 || psi.length != 1)
                throw new RuntimeException("Birth, death and sampling rate have to be 1-dimensional!");
        }

        T = origin.get().getValue();
        ntaxa = m_tree.get().getLeafNodeCount();

    }


    @Override
    public Double updateRatesAndTimes(Tree tree){

        super.updateRatesAndTimes(tree);

        T = origin.get().getValue();
        ntaxa = tree.getLeafNodeCount();

        S0 = (S0_input.get().getArrayValue());

        dS = m_dS.get().getValues();

        dE = (m_dE.get() !=null)? m_dE.get().getValues(): (new Double[dS.length]);
        if (dE[0]==null) Arrays.fill(dE,0.);

        dR = m_dR.get().getValues();


        double cumS = S0 - 1 ;
        double b = birth[0]/S0 ;
        double time;

        Double[] birthSIR = new Double[dim];
        double I = 1.;
        double R = 0.;

        birthSIR[0] = b * cumS;
        for (int i = 0; i < dim-1; i++){

            cumS -= dS[i];
            birthSIR[i+1] = b * cumS;

            I += dS[i] - dE[i] - dR[i];
            R += dR[i];
            time = (i+1.)/dim * T;

            if ( checkTreeConsistent.get() && (I<=0. || I < lineageCountAtTime(T-time, tree)))
                return Double.NEGATIVE_INFINITY;

        }

        if (cumS < 0 || S0 - cumS < m_tree.get().getLeafNodeCount() || S0 != (cumS+I+R))
            return Double.NEGATIVE_INFINITY;

        adjustBirthRates(birthSIR);
        return 0.;

    }


    @Override
    public Boolean isBDSIR(){
        return true;
    }

    public int getSIRdimension(){
        return dim;
    }
}
