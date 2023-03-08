import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import string
import pyphi

CAUSE = pyphi.direction.Direction.CAUSE  
EFFECT = pyphi.direction.Direction.EFFECT 

#### EXAMPLE SYSTEMS

def mexican_hat_grid(
    n_nodes=4, self_loop=0.95, close=0.1, far=-0.025, threshold=1.0, determinism=10, periodic=True, plot_activation_function=True
):

    # local functions
    # Create weight matrix (nearest neighbour)
    def get_weight_matrix(n_nodes, periodic):
        weight = [close, self_loop, close]
        weights = np.ones((n_nodes, n_nodes))*far
        for n in range(n_nodes):
            weights[n, n - 1] = weight[0]
            weights[n, n] = weight[1]
            weights[n, int((n + 1) % n_nodes)] = weight[2] if periodic else 0

        return weights

    # Logistic function
    def LogFunc(x, determinism, threshold):
        y = 1 / (1 + np.e ** (-determinism * (x - threshold)))
        return y

    # produce weight matrix
    weights = get_weight_matrix(n_nodes, periodic)
        
    # producing transition probability matrix
    tpm = [
            [
                LogFunc(
                    sum(state * np.array([weights[z][n] for n in range(n_nodes)])),
                    determinism,
                    threshold,
                )
                for z in range(n_nodes)
            ]
            for state in pyphi.utils.all_states(n_nodes)
        ]
    
    if plot_activation_function:
        
        input_state_weights = [
            sum(state * np.array([weights[1][n] for n in range(3)]))
            for state in pyphi.utils.all_states(3)
        ]
        input_state_P = np.array([LogFunc(xx, determinism, threshold) for xx in input_state_weights])
        
        x = np.linspace(min(input_state_weights)-1, max(input_state_weights)+1, 1001)
        y = np.array([LogFunc(xx, determinism, threshold) for xx in x])
        
        plt.plot(x,y,'k')
        plt.scatter(input_state_weights, input_state_P)
        
        plt.title('Activation Function for units\n(dots show realizable input states)')
        plt.xlabel('Summed input weight to unit')
        plt.ylabel('Probability of activating')
        
    return pyphi.Network(
        tpm=tpm,
        cm=np.array(
            [[float(1) if w else 0 for w in weights[n]] for n in range(len(weights))]
        ).T,
        node_labels=[string.ascii_uppercase[n] for n in range(len(weights))],
    )

def grid_1d_nearest_neighbor(
    n_nodes=4, self_loop=0.95, lateral=0.025, threshold=1.0, determinism=10, periodic=True, plot_activation_function=True
):

    # local functions
    # Create weight matrix (nearest neighbour)
    def get_weight_matrix(n_nodes, periodic):
        weight = [lateral, self_loop, lateral]
        weights = np.zeros((n_nodes, n_nodes))
        for n in range(n_nodes):
            weights[n, n - 1] = weight[0]
            weights[n, n] = weight[1]
            weights[n, int((n + 1) % n_nodes)] = weight[2] if periodic else 0

        return weights

    # Logistic function
    def LogFunc(x, determinism, threshold):
        y = 1 / (1 + np.e ** (-determinism * (x - threshold)))
        return y

    # produce weight matrix
    weights = get_weight_matrix(n_nodes, periodic)
        
    # producing transition probability matrix
    tpm = [
            [
                LogFunc(
                    sum(state * np.array([weights[z][n] for n in range(n_nodes)])),
                    determinism,
                    threshold,
                )
                for z in range(n_nodes)
            ]
            for state in pyphi.utils.all_states(n_nodes)
        ]
    
    if plot_activation_function:
        
        input_state_weights = [
            sum(state * np.array([weights[1][n] for n in range(3)]))
            for state in pyphi.utils.all_states(3)
        ]
        input_state_P = np.array([LogFunc(xx, determinism, threshold) for xx in input_state_weights])
        
        x = np.linspace(min(input_state_weights)-1, max(input_state_weights)+1, 1001)
        y = np.array([LogFunc(xx, determinism, threshold) for xx in x])
        
        plt.plot(x,y,'k')
        plt.scatter(input_state_weights, input_state_P)
        
        plt.title('Activation Function for units\n(dots show realizable input states)')
        plt.xlabel('Summed input weight to unit')
        plt.ylabel('Probability of activating')
        
    return pyphi.Network(
        tpm=tpm,
        cm=np.array(
            [[float(1) if w else 0 for w in weights[n]] for n in range(len(weights))]
        ).T,
        node_labels=[string.ascii_uppercase[n] for n in range(len(weights))],
    )


#### UTILS 
def get_unit_tpm(network, unit_id):

    # run "experiment"
    effect_of_purview_state = dict()

    for state in pyphi.utils.all_states(len(network.node_indices)):
        purview_state = tuple(
            [s for i, s in enumerate(state) if i==unit_id]
        )
        effect = network.tpm[state]

        if purview_state in effect_of_purview_state:
            effect_of_purview_state[purview_state].append(list(effect))
        else:
            effect_of_purview_state[purview_state] = [list(effect)]

    # construct tpm
    effect_tpm = np.array(
        [
            list(np.array(effect_of_purview_state[purview_state]).mean(0))
            for purview_state in pyphi.utils.all_states(1)
        ]
    )

    em = np.squeeze(effect_tpm)[:, unit_id][..., np.newaxis]
    sbs_effect_tpm = np.concatenate((1 - em, em), 1)
    
    labels = [
        l for i, l in enumerate(network.node_labels) if i==unit_id
    ]
    
    # get tpm df
    df_tpm = show_tpm(
        sbs_effect_tpm,
        labels,
        labels,
        style="pandas",
        caption=f'Unit TPM for {labels}.',
    )

    return df_tpm


def get_mechanism_effect_tpm(subsystem, mechanism):

    # get purview
    possible_output_elements = tuple(
        set(
            [
                elem
                for purview in subsystem.potential_purviews(
                    EFFECT, mechanism
                )
                for elem in purview
            ]
        )
    )

    # run "experiment"
    effect_of_mechanism_state = dict()

    for state in pyphi.utils.all_states(len(subsystem.node_indices)):
        mechanism_state = tuple([s for i, s in enumerate(state) if i in mechanism])
        effect = subsystem.tpm[state]

        if mechanism_state in effect_of_mechanism_state:
            effect_of_mechanism_state[mechanism_state].append(list(effect))
        else:
            effect_of_mechanism_state[mechanism_state] = [list(effect)]

    # construct tpm
    effect_tpm = np.array(
        [
            list(np.array(effect_of_mechanism_state[mechanism_state]).mean(0))
            for mechanism_state in pyphi.utils.all_states(len(mechanism))
        ]
    )

    sbs_effect_tpm = pyphi.convert.sbn2sbs(
        np.squeeze(effect_tpm)[:, possible_output_elements]
    )

    past_labels = [l for i, l in enumerate(subsystem.node_labels) if i in mechanism]
    future_labels = [
        l for i, l in enumerate(subsystem.node_labels) if i in possible_output_elements
    ]
    # get tpm df
    df_tpm = show_tpm(
        sbs_effect_tpm,
        past_labels,
        future_labels,
        style="pandas",
        caption=f'Effect TPM from the present state of {"".join(past_labels)} to the future state of {"".join(future_labels)}.',
    )

    return df_tpm


def get_mechanism_cause_tpm(subsystem, mechanism):

    # get purview
    possible_input_elements = tuple(
        set(
            [
                elem
                for purview in subsystem.potential_purviews(
                    CAUSE, mechanism
                )
                for elem in purview
            ]
        )
    )

    # run "experiment"
    effect_of_purview_state = dict()

    for state in pyphi.utils.all_states(len(subsystem.node_indices)):
        purview_state = tuple(
            [s for i, s in enumerate(state) if i in possible_input_elements]
        )
        effect = subsystem.tpm[state]

        if purview_state in effect_of_purview_state:
            effect_of_purview_state[purview_state].append(list(effect))
        else:
            effect_of_purview_state[purview_state] = [list(effect)]

    # construct tpm
    effect_tpm = np.array(
        [
            list(np.array(effect_of_purview_state[purview_state]).mean(0))
            for purview_state in pyphi.utils.all_states(len(possible_input_elements))
        ]
    )

    if len(mechanism) == 1:
        em = np.squeeze(effect_tpm)[:, mechanism]
        sbs_effect_tpm = np.concatenate((1 - em, em), 1)
    else:
        sbs_effect_tpm = pyphi.convert.sbn2sbs(np.squeeze(effect_tpm)[:, mechanism])
    past_labels = [
        l for i, l in enumerate(subsystem.node_labels) if i in possible_input_elements
    ]
    future_labels = [l for i, l in enumerate(subsystem.node_labels) if i in mechanism]
    # get tpm df
    df_tpm = show_tpm(
        sbs_effect_tpm,
        past_labels,
        future_labels,
        style="pandas",
        caption=f'Cause TPM from the past state of {"".join(past_labels)} to the present state of {"".join(future_labels)}.',
    )

    return df_tpm



def get_mechanism_tpm(subsystem, mechanism, direction=CAUSE):
    if direction == CAUSE:
        return get_mechanism_cause_tpm(subsystem, mechanism)
    else:
        return get_mechanism_effect_tpm(subsystem, mechanism)

    
def get_candidate_cause_or_effect(subsystem, mechanism, purview, direction, state, partition):
    phi, partitioned_repertoire = subsystem.evaluate_partition(
            direction,
            mechanism,
            purview,
            partition,
            state=state,
    )
    
    ria = pyphi.models.mechanism.RepertoireIrreducibilityAnalysis(
        phi,
        direction,
        mechanism,
        purview,
        partition,
        subsystem.repertoire(
            direction, mechanism, purview
        ),
        partitioned_repertoire,
        specified_index=None,
        specified_state=np.array([state]),
        mechanism_state=None,
        purview_state=None,
        node_labels=subsystem.node_labels,
    )
    ria.ria = ria
    return pyphi.models.MaximallyIrreducibleCauseOrEffect(ria)



def get_candidate_distinction(
    subsystem,
    mechanism,
    candidate_cause_purview=False,
    candidate_effect_purview=False,
    candidate_cause_state=False,
    candidate_effect_state=False,
    candidate_cause_cut=False,
    candidate_effect_cut=False,
):

    # setting default args, if not given
    
    # the default purviews are the maximally irreducible ones
    if not candidate_cause_purview:
        mic = subsystem.find_mice(CAUSE, mechanism)
        candidate_cause_purview = mic.purview
    if not candidate_effect_purview:
        mie = subsystem.find_mice(EFFECT, mechanism)
        candidate_effect_purview = mie.purview
        
    if candidate_cause_purview==() or candidate_effect_purview==():
        return 
    
    # default states are the maximally informative ones
    if not candidate_cause_state:
        candidate_cause_state = subsystem.find_maximal_state_under_complete_partition(CAUSE, mechanism, candidate_cause_purview)[0]
    if not candidate_effect_state:
        candidate_effect_state = subsystem.find_maximal_state_under_complete_partition(EFFECT, mechanism, candidate_effect_purview)[0]

    # the default cut is the minimum information partition
    if not candidate_cause_cut:
        if mic in locals():
            candidate_cause_cut = mic.mip
        else:
            candidate_cause_cut = subsystem.find_mice(CAUSE, mechanism).mip
            
    if not candidate_effect_cut:
        if mie in locals():
            candidate_effect_cut = mie.mip
        else:
            candidate_effect_cut = subsystem.find_mice(EFFECT, mechanism).mip
        
    # Computing the cause and the effect
    candidate_cause = get_candidate_cause_or_effect(
        subsystem,
        mechanism,
        candidate_cause_purview,
        CAUSE,
        candidate_cause_state,
        candidate_cause_cut,
    )
    candidate_effect = get_candidate_cause_or_effect(
        subsystem,
        mechanism,
        candidate_effect_purview,
        EFFECT,
        candidate_effect_state,
        candidate_effect_cut,
    )

    candidate_distinction = pyphi.models.Concept(
        mechanism=mechanism, cause=candidate_cause, effect=candidate_effect, subsystem=subsystem
    )
    return candidate_distinction



##### Visualization

def show_tpm(tpm, past_labels=None, future_labels=None, style="pandas", caption=None):
    
    if len(tpm.shape)>2:
        tpm = pyphi.convert.sbn2sbs(tpm)

    
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    future_states = list(pyphi.utils.all_states(len(future_labels)))

    assert (
        len(past_states),
        len(future_states),
    ) == tpm.shape, "mismatch in TPM shape and labels provided."

    if style == "pandas":

        columns = pd.MultiIndex.from_tuples(future_states, names=future_labels)
        index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

        TPM = pd.DataFrame(tpm, index=index, columns=columns)
        (TPM.style.background_gradient(cmap="gray_r", low=0, high=1, axis=None)
            .format(precision=2)
            .set_caption(
                f'TPM from the past state of {"".join(past_labels)} to the future state of {"".join(future_labels)}.'
                if caption == None
                else caption
            ))
        return TPM

def highlight_cell(col, col_label, row_label, color="lightblue"):
   # check if col is a column we want to highlight
    if col.name == col_label:
        # a boolean mask where True represents a row we want to highlight
        mask = (col.index == row_label)
        # return a list of string styles (e.g. ["", "background-color: yellow"])
        return [
            'background-color: {}'.format(color)
            if val_bool else ""
            for val_bool in mask
        ]
    else:
        # return an array of empty strings that has the same size as col (e.g. ["",""])
        return np.full_like(col, "", dtype="str")

def highlight_transition_probability(TPM,input_state, output_state, color="lightblue"):
    return TPM.style.apply(highlight_cell, col_label=output_state, row_label=input_state, color=color)
    
def network_tpm(network):
    return show_tpm(
        network.tpm, 
        past_labels=network.node_labels, 
        future_labels=network.node_labels, 
        style="pandas", 
        caption='The full TPM of the substrate {}'.format([l for l in network.node_labels])
    )
    
def plot_repertoire(ax, repertoire, states, stagger=0, color="black"):

    xaxis = np.array(range(len(repertoire)))

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax.spines["bottom"].set_position("zero")

    # Eliminate upper and right axes
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")

    ax.bar(xaxis - stagger, repertoire, width=0.7, color=color)

    ax.set_xticks(xaxis, purview_states, rotation=90)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylim([-0.1, 1.1])

    return ax


def cause_repertoire(subsystem, mechanism, purview):

    repertoire = pyphi.distribution.flatten(
        np.squeeze(subsystem.cause_repertoire(mechanism, purview))
    )
    purview_states = pyphi.utils.all_states(len(purview))

    # plotting first histogram
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = plot_repertoire(ax, repertoire, purview_states, stagger=0.0)
    return fig


def effect_repertoire(subsystem, mechanism, purview):

    repertoire = pyphi.distribution.flatten(
        np.squeeze(subsystem.effect_repertoire(mechanism, purview))
    )
    purview_states = pyphi.utils.all_states(len(purview))

    # plotting first histogram
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = plot_repertoire(ax, repertoire, purview_states, stagger=0.0)
    return fig


def cause_repertoire_integration(subsystem, mechanism, purview):

    repertoire = pyphi.distribution.flatten(
        np.squeeze(subsystem.cause_repertoire(mechanism, purview))
    )

    direction = CAUSE
    cut_repertoire = pyphi.distribution.flatten(
        np.squeeze(
            subsystem.find_mip(
                direction,
                mechanism,
                purview,
                state=subsystem.find_maximal_state_under_complete_partition(
                    direction, mechanism, purview
                )[0],
            ).partitioned_repertoire
        )
    )
    purview_states = pyphi.utils.all_states(len(purview))

    # plotting first histogram
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = plot_repertoire(ax, cut_repertoire, purview_states, stagger=0.1, color="red")
    ax = plot_repertoire(ax, repertoire, purview_states, stagger=0.0)
    return fig


def effect_repertoire_integration(subsystem, mechanism, purview):

    repertoire = pyphi.distribution.flatten(
        np.squeeze(subsystem.cause_repertoire(mechanism, purview))
    )

    direction = EFFECT
    cut_repertoire = pyphi.distribution.flatten(
        np.squeeze(
            subsystem.find_mip(
                direction,
                mechanism,
                purview,
                state=subsystem.find_maximal_state_under_complete_partition(
                    direction, mechanism, purview
                )[0],
            ).partitioned_repertoire
        )
    )
    purview_states = pyphi.utils.all_states(len(purview))

    # plotting first histogram
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = plot_repertoire(ax, cut_repertoire, purview_states, stagger=0.1, color="green")
    ax = plot_repertoire(ax, repertoire, purview_states, stagger=0.0)
    return fig


def cause_repertoire_information(subsystem, mechanism, purview):

    repertoire = pyphi.distribution.flatten(
        np.squeeze(subsystem.cause_repertoire(mechanism, purview))
    )

    direction = CAUSE
    cut_repertoire = pyphi.distribution.flatten(
        np.squeeze(
            subsystem.unconstrained_repertoire(
                direction,
                purview,
                ),
            )
        )
    purview_states = pyphi.utils.all_states(len(purview))

    # plotting first histogram
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = plot_repertoire(ax, cut_repertoire, purview_states, stagger=0.1, color="gray")
    ax = plot_repertoire(ax, repertoire, purview_states, stagger=0.0)
    return fig


def effect_repertoire_information(subsystem, mechanism, purview):

    repertoire = pyphi.distribution.flatten(
        np.squeeze(subsystem.cause_repertoire(mechanism, purview))
    )

    direction = EFFECT
    cut_repertoire = pyphi.distribution.flatten(
        np.squeeze(
            subsystem.unconstrained_repertoire(
                direction,
                purview,
                ),
            )
        )
    purview_states = pyphi.utils.all_states(len(purview))

    # plotting first histogram
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax = plot_repertoire(ax, cut_repertoire, purview_states, stagger=0.1, color="gray")
    ax = plot_repertoire(ax, repertoire, purview_states, stagger=0.0)
    return fig