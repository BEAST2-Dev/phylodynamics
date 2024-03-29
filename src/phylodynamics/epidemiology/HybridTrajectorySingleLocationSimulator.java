package phylodynamics.epidemiology;

import beast.base.inference.Operator;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Function;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;


/**
 * User: Denise
 * Date: Jun 29, 2012
 * Time: 7:19:46 PM
 */
@Description("Operator that simulates new S(E)IR trajectory")
public class HybridTrajectorySingleLocationSimulator extends Operator {


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
    public Input<RealParameter> origin =
            new Input<RealParameter>("origin", "The origin of infection (> treeheight)", Input.Validate.REQUIRED);

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


    @Override
    public void initAndValidate() {

        S0 = (int) (S0_input.get().getArrayValue());

        if (birth.get() != null && death.get() != null && sampling.get() != null){

            b = birth.get().getArrayValue();

            int dim = death.get().getDimension();
            if (dim != sampling.get().getDimension()) throw new RuntimeException("Error: Death and sampling must have equalt dimensions!");
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

        ntaxa = tree.getLeafNodeCount();

        // initialize trajectory
        current =  SIR.get();
        if ( !current.initTraj(new Double[]{b}, 0., d, s, 0., T, ntaxa, 100, m, current.times, false))
            throw new RuntimeException("Could not find suitable trajectory. Please try different epi parameters!");

        dS_input.get().assignFromWithoutID(new RealParameter(current.dS));
        dR_input.get().assignFromWithoutID(new RealParameter(current.dR));
    }




    @Override
    public double proposal() {

        S0 = (int) (S0_input.get().getArrayValue());

        b = birth.get().getArrayValue();

        int dim = death.get().getDimension();
        if (dim != sampling.get().getDimension()) throw new RuntimeException("Error: Death and sampling must have equalt dimensions!");
        d = new Double[dim];
        s = new Double[dim];

        for (int i = 0; i<dim; i++){
            d[i] = death.get().getArrayValue(i);
            s[i] = sampling.get().getArrayValue(i);
        }

        T =  origin.get().getValue();

        // simulate trajectory
        current =  SIR.get();
        if ( !current.refresh(S0,0, new Double[]{b}, 0., d, s, 0., T, ntaxa, 1, m, current.times, false) /* && Randomizer.nextDouble()>.75*/){
            return Double.NEGATIVE_INFINITY;
        }

        for (int i =0; i < m-1; i++){
            dS_input.get().setValue(i, current.dS[i]);
            dR_input.get().setValue(i, current.dR[i]);
        }

        return 0.;
    }

    public static void main(String[] args){
    }
}
