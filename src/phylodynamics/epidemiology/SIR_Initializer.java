package phylodynamics.epidemiology;


import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import beast.base.inference.util.RPNcalculator;

import java.util.List;


/**
 * @author dkuh004
 *         Date: Apr 12, 2013
 *         Time: 10:46:37 AM
 */
@Description("Initializes a stichastic SIR trajectory")
public class SIR_Initializer extends BEASTObject implements StateNodeInitialiser {



    public Input<HybridSEIREpidemic> SIR =
            new Input<HybridSEIREpidemic>("SIR", "SIR trajectory calculation node - simulates dS", Input.Validate.REQUIRED);

    public Input<Function> S0_input =
            new Input<Function>("S0", "The numbers of susceptible individuals", Input.Validate.REQUIRED);

    public Input<RealParameter> dS_input =
            new Input<RealParameter>("dS", "dS vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);
    public Input<RealParameter> dR_input =
            new Input<RealParameter>("dR", "dR vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);

    public Input<Function> birth =
            new Input<Function>("birth", "birth rate vector with rate per location",  Input.Validate.REQUIRED);
    public Input<Function> death =
            new Input<Function>("death", "death rate vector with rate per location",  Input.Validate.REQUIRED);
    public Input<Function> sampling =
            new Input<Function>("sampling", "sampling rate vector with rate per location",  Input.Validate.REQUIRED);

    public Input<Tree> m_tree =
            new Input<Tree>("tree", "The phylogenetic tree being estimated",  Input.Validate.REQUIRED);
    public Input<RealParameter> orig_root =
            new Input<RealParameter>("orig_root", "The origin of infection x0", Input.Validate.REQUIRED);

    public Input<RealParameter> R0 =
            new Input<RealParameter>("R0", "The basic reproduction number");
    public Input<RealParameter> becomeUninfectiousRate =
            new Input<RealParameter>("becomeUninfectiousRate", "Rate at which individuals become uninfectious (throuch recovery or sampling)");
    public Input<RealParameter> samplingProportion =
            new Input<RealParameter>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate");



    int m;
    Integer S0;
    int ntaxa;

    Tree tree;
    double T;
    HybridSEIREpidemic current;

    Double b;
    Double[] d;
    Double[] s;
    Boolean birthChanges;
    Boolean deathChanges;
    Boolean samplingChanges;

    double scaler;


    @Override
    public void initAndValidate(){
         scaler = 1.;
    }

    @Override
    public void initStateNodes() {


        S0 = (int) (S0_input.get().getArrayValue());

        if (birth.get() != null && death.get() != null && sampling.get() != null){

            b = birth.get().getArrayValue() *(1.+scaler);

            int dim = death.get().getDimension();
            if (dim != sampling.get().getDimension()) throw new RuntimeException("Error: Death and sampling must have equalt dimensions!");
            d = new Double[dim];
            s = new Double[dim];

            for (int i = 0; i<dim; i++){
                d[i] = death.get().getArrayValue(i)  *(1.-scaler);
                s[i] = sampling.get().getArrayValue(i) ;
            }
        }

        else{
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }

        if (birth.get() instanceof RPNcalculator && (R0.get() == null || becomeUninfectiousRate.get()==null || samplingProportion.get()==null))
            throw new RuntimeException("R0, becomeUninfectiousRate and samplingProportion need to be specified");



        m = SIR.get().Nsamples.get();

        tree = m_tree.get();
        T = tree.getRoot().getHeight() + orig_root.get().getValue();

        ntaxa = tree.getLeafNodeCount();

        // initialize trajectory
        current =  SIR.get();
//        if ( !current.init(b, 0., d, s, T, ntaxa, 100, m, current.times))
//            throw new RuntimeException("Could not find suitable trajectory. Please try different epi parameters!");

        if (!current.initTraj(new Double[]{b}, 0., d, s, 0., T, ntaxa, 100, m, current.times, false)){

            int count = 0;

            do  {

                scaler = Randomizer.nextDouble() ; // if first initialization didn't work, try tweaking the rates a little
                b*= (1.+scaler);
                d[0]*= (1.-scaler);
                count++;

            } while (count<100 && !current.initTraj(new Double[]{b}, 0., d, s, 0., T, ntaxa, 100, m, current.times, false));
        }

        dS_input.get().assignFromWithoutID(new RealParameter(current.dS));
        dR_input.get().assignFromWithoutID(new RealParameter(current.dR));

        if (R0.get()!=null) R0.get().setValue(0,b/(d[0]+s[0]));
        else ((RealParameter) birth.get()).setValue(0,b);

        if (becomeUninfectiousRate.get()!=null) becomeUninfectiousRate.get().setValue(0,(d[0]+s[0]));
        else ((RealParameter) death.get()).setValue(0,d[0]);

        if (samplingProportion.get()!=null) samplingProportion.get().setValue(0,s[0]/(d[0]+s[0]));
        else ((RealParameter) sampling.get()).setValue(0,s[0]);


        scaler = Randomizer.nextDouble() ; // if first initialization didn't work, try tweaking the rates a little

    }


    /**
     * @return list of StateNodes that are initialised
     */
    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {

        stateNodes.add(dS_input.get());
        stateNodes.add(dR_input.get());
    }

}
