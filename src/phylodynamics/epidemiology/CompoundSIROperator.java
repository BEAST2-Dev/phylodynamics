package phylodynamics.epidemiology;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.CalculationNode;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;


/**
 * User: Denise
 * Date: Jan 4, 2013
 * Time: 12:17:34 PM
 */
@Description("Operator that simulates new S(E)IR trajectory")
public class CompoundSIROperator extends Operator {

    public Input<Operator> affectingOperator =
            new Input<Operator>("affectingOperator", "Operator that affects S(E)IR trajectory", Input.Validate.REQUIRED);

    public Input<HybridSEIREpidemic> SIR =
            new Input<HybridSEIREpidemic>("SIR", "SIR trajectory calculation node - simulates dS", Input.Validate.REQUIRED);

    public Input<Function> S0_input =
            new Input<Function>("S0", "The numbers of susceptible individuals", Input.Validate.REQUIRED);

    public Input<RealParameter> dS_input =
            new Input<RealParameter>("dS", "dS vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);
    public Input<RealParameter> dR_input =
            new Input<RealParameter>("dR", "dR vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);

    public Input<RealParameter> susceptiblePopulationInput = new Input<RealParameter>("susceptiblePopulation", "susceptible population size at times eventTimes");

    public Input<RealParameter> infectedPopulationInput =
            new Input<RealParameter>("infectedPopulation", "infectedPopulation over time (relative to eventTimes)", Input.Validate.XOR, dS_input);

    public Input<RealParameter> eventTimeInput =
            new Input<RealParameter>("eventTimes", "ordered times at which events in SIR happen", Input.Validate.XOR, dR_input);

    public Input<Function> birth =
            new Input<Function>("birth", "birth rate vector with rate per location",  Input.Validate.REQUIRED);
    public Input<Function> death =
            new Input<Function>("death", "death rate vector with rate per location",  Input.Validate.REQUIRED);
    public Input<Function> sampling =
            new Input<Function>("sampling", "sampling rate vector with rate per location",  Input.Validate.REQUIRED);

    public Input<RealParameter> loseImmunityRate = new Input<RealParameter>("loseImmunityRate", "The rate at which recovereds lose immunity");

    public Input<Tree> m_tree =
            new Input<Tree>("tree", "The phylogenetic tree being estimated",  Input.Validate.REQUIRED);
    public Input<RealParameter> origin =
            new Input<RealParameter>("origin", "The origin of infection x0", Input.Validate.REQUIRED);

    int m;
    Integer S0;
    int ntaxa;

    Tree tree;
    double T;
    HybridSEIREpidemic current;

    Double[] b;
    Double[] d;
    Double[] s;
    Boolean birthChanges;
    Boolean deathChanges;
    Boolean samplingChanges;


    @Override
    public void initAndValidate() {


        affectingOperator.get().initAndValidate();

        S0 = (int) (S0_input.get().getArrayValue());

        if (birth.get() != null && death.get() != null && sampling.get() != null){

            b = new Double[birth.get().getDimension()];

            for (int i=0; i<b.length; i++){
                b[i] = birth.get().getArrayValue(i);
            }


            int dim = death.get().getDimension();
            if (dim != sampling.get().getDimension()) throw new RuntimeException("Error: Death and sampling must have equal dimensions!");
            d = new Double[dim];
            s = new Double[dim];

            for (int i = 0; i<dim; i++){
                d[i] = death.get().getArrayValue(i);
                s[i] = sampling.get().getArrayValue(i);
            }
        }

        else{
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }


        m = SIR.get().Nsamples.get();

        tree = m_tree.get();
        T = origin.get().getValue();
        if (T<=0) T=tree.getRoot().getHeight();

        ntaxa = tree.getLeafNodeCount();

        // initialize trajectory
        current =  SIR.get();
        if ( !current.initTraj(b, 0., d, s, (loseImmunityRate.get()==null)?0.:loseImmunityRate.get().getValue(), T, ntaxa, 100, m, current.times, dS_input.get()==null))
            throw new RuntimeException("Could not find suitable trajectory. Please try different epi parameters!");

        if (dS_input.get() !=null){

            dS_input.get().assignFromWithoutID(new RealParameter(current.dS));
            dR_input.get().assignFromWithoutID(new RealParameter(current.dR));
        }
        else {
            susceptiblePopulationInput.get().assignFromWithoutID(new RealParameter(current.S));
            infectedPopulationInput.get().assignFromWithoutID(new RealParameter(current.I));
            eventTimeInput.get().assignFromWithoutID(new RealParameter(current.eventTimes));

        }
    }

    @Override
    public double proposal() {

        double hastingsRatio = affectingOperator.get().proposal();

        S0 = (int) (S0_input.get().getArrayValue());

//        b = new Double[birth.get().getDimension()];

        for (int i=0; i<b.length; i++){
            b[i] = birth.get().getArrayValue(i);
        }

        int dim = death.get().getDimension();
        if (dim != sampling.get().getDimension()) throw new RuntimeException("Error: Death and sampling must have equal dimensions!");
//        d = new Double[dim];
//        s = new Double[dim];

        for (int i = 0; i<dim; i++){
            d[i] = death.get().getArrayValue(i);
            s[i] = sampling.get().getArrayValue(i);
        }

        T = origin.get().getValue();
        if (T<=0) T=tree.getRoot().getHeight();

        // simulate trajectory
        current =  SIR.get();
        if ( !current.refresh(S0,0, b, 0., d, s,(loseImmunityRate.get()==null)?0.:loseImmunityRate.get().getValue(), T, ntaxa, 1, m, current.times, dS_input.get()==null) /* && Randomizer.nextDouble()>.75*/){
            return Double.NEGATIVE_INFINITY;
        }

        if (dS_input.get()!=null){
            for (int i =0; i < m-1; i++){
                dS_input.get().setValue(i, current.dS[i]);
                dR_input.get().setValue(i, current.dR[i]);
            }
        }
        else{
            for (int i =0; i < current.eventTimes.length; i++){
                susceptiblePopulationInput.get().setValue(i, current.S[i]);
                infectedPopulationInput.get().setValue(i, current.I[i]);
                eventTimeInput.get().setValue(i, current.eventTimes[i]);
            }

        }
        return hastingsRatio;
    }


    /**
     * return list of state nodes that this operator operates on.
     * state nodes that are input to the operator but are never changed
     * in a proposal should not be listed
     */
    @Override
    public List<StateNode> listStateNodes() {

        final List<StateNode> list = affectingOperator.get().listStateNodes();

        list.add(dS_input.get());
        list.add(dR_input.get());

        return list;
    }

    public static void main(String[] args){
    }
}
