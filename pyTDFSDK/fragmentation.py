import pandas as pd


def has_fragmentation_tables(analysis):
    """
    Check whether the database has the Fragmentation* tables of schema v5 and later. In schema v5 the older MS/MS
    tables (i.e. PasefFrameMsMsInfo, FrameMsMsInfo) are views that join the CID settings and so drop frames that use
    other fragmentation (i.e. ExD). The Fragmentation* tables have to be read instead.

    :param analysis: Dictionary of tables from pyTDFSDK.classes.TdfData.analysis or
        pyTDFSDK.classes.TsfData.analysis.
    :type analysis: dict
    :return: True if the Fragmentation* tables exist.
    :rtype: bool
    """
    return 'FragmentationStep' in analysis


def get_fragmentation_steps(analysis):
    """
    Get one row per fragmentation step with its isolation window and activation settings. The CID step that follows
    an electron based step (ExD then CID, i.e. ECcID) is linked through FragmentationStepToStep; its collision energy
    is used when the step has no collision energy of its own.

    :param analysis: Dictionary of tables from pyTDFSDK.classes.TdfData.analysis or
        pyTDFSDK.classes.TsfData.analysis.
    :type analysis: dict
    :return: Table with the columns StepId, Precursor, IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd,
        CollisionEnergy, ReactionTime, ElectronEnergy, Parent, PrecursorCharge and PrecursorGroup.
    :rtype: pandas.DataFrame
    """
    steps = analysis['FragmentationStep'].rename(columns={'Id': 'StepId'})
    isolation = analysis['FragmentationIsolationWindow'].rename(columns={'Id': 'IsolationWindow'})
    steps = steps.merge(isolation, on='IsolationWindow', how='left')

    cid = analysis['FragmentationCid'].rename(columns={'FragmentationStepId': 'StepId'})
    exd = analysis['FragmentationExd'].rename(columns={'FragmentationStepId': 'StepId'})
    etd = analysis['FragmentationEtd'].rename(columns={'FragmentationStepId': 'StepId',
                                                       'ReactionTime': 'EtdReactionTime'})
    steps = steps.merge(cid, on='StepId', how='left').merge(exd, on='StepId', how='left')
    steps = steps.merge(etd, on='StepId', how='left')

    link = analysis['FragmentationStepToStep'].rename(columns={'FragmentationStepId': 'StepId',
                                                               'FragmentationStepOriginId': 'OriginId'})
    linked_cid = link.merge(cid, on='StepId')[['OriginId', 'CollisionEnergy']]
    linked_cid = linked_cid.rename(columns={'OriginId': 'StepId', 'CollisionEnergy': 'LinkedCollisionEnergy'})
    steps = steps.merge(linked_cid.drop_duplicates('StepId'), on='StepId', how='left')
    steps['CollisionEnergy'] = steps['CollisionEnergy'].fillna(steps['LinkedCollisionEnergy'])
    steps['ReactionTime'] = steps['ReactionTime'].fillna(steps['EtdReactionTime'])

    # PrecursorGroup was added in schema 5.1 (guided PASEF: links the CID and ExD precursors of one ion).
    precursors = analysis['FragmentationPrecursor']
    columns = ['Id', 'Parent', 'Charge'] + (['PrecursorGroup'] if 'PrecursorGroup' in precursors.columns else [])
    precursors = precursors[columns].rename(columns={'Id': 'Precursor', 'Charge': 'PrecursorCharge'})
    steps = steps.merge(precursors, on='Precursor', how='left')
    if 'PrecursorGroup' not in steps.columns:
        steps['PrecursorGroup'] = None

    return steps[['StepId', 'Precursor', 'IsolationMz', 'IsolationWidth', 'ScanNumBegin', 'ScanNumEnd',
                  'CollisionEnergy', 'ReactionTime', 'ElectronEnergy', 'Parent', 'PrecursorCharge',
                  'PrecursorGroup']]


def get_frame_steps(analysis):
    """
    Get one row per frame and fragmentation step, ordered by step ID.

    :param analysis: Dictionary of tables from pyTDFSDK.classes.TdfData.analysis or
        pyTDFSDK.classes.TsfData.analysis.
    :type analysis: dict
    :return: Table with the column Frame and the columns of pyTDFSDK.fragmentation.get_fragmentation_steps().
    :rtype: pandas.DataFrame
    """
    experiments = analysis['FragmentationExperiment'].rename(columns={'FragmentationStep': 'StepId'})
    frame_steps = analysis['FragmentationFrameToExperiment'].merge(experiments, on='GroupId')
    frame_steps = frame_steps.merge(get_fragmentation_steps(analysis), on='StepId', how='left')
    return frame_steps.sort_values('StepId', kind='stable').reset_index(drop=True)


def _group_records(frame_steps, key):
    """
    Group rows of pyTDFSDK.fragmentation.get_frame_steps() by a column. Missing values are None.
    """
    frame_steps = frame_steps[frame_steps[key].notna()]
    frame_steps = frame_steps.astype(object).where(frame_steps.notna(), None)
    return {int(value): rows.to_dict(orient='records') for value, rows in frame_steps.groupby(key)}


def get_precursor_steps(data):
    """
    Get the fragmentation steps of each precursor as a dictionary of precursor ID to a list of rows (one row for each
    frame the precursor was fragmented in, ordered by step ID). The first row holds the settings of the precursor.
    The result is stored in the data object so it is only built once.

    :param data: Object with the analysis tables.
    :type data: pyTDFSDK.classes.TdfData | pyTDFSDK.classes.TsfData
    :return: Dictionary of precursor ID to list of dictionaries.
    :rtype: dict
    """
    if getattr(data, 'precursor_steps', None) is None:
        data.precursor_steps = _group_records(get_frame_steps(data.analysis), 'Precursor')
    return data.precursor_steps


def get_steps_by_frame(data):
    """
    Get the fragmentation steps of each frame as a dictionary of frame ID to a list of rows ordered by step ID. The
    result is stored in the data object so it is only built once.

    :param data: Object with the analysis tables.
    :type data: pyTDFSDK.classes.TdfData | pyTDFSDK.classes.TsfData
    :return: Dictionary of frame ID to list of dictionaries.
    :rtype: dict
    """
    if getattr(data, 'steps_by_frame', None) is None:
        data.steps_by_frame = _group_records(get_frame_steps(data.analysis), 'Frame')
    return data.steps_by_frame


def get_activation(step):
    """
    Get the activation method of a fragmentation step.

    :param step: Row from pyTDFSDK.fragmentation.get_precursor_steps() or get_steps_by_frame(), or a row of an
        older MS/MS table (which has no ElectronEnergy or ReactionTime).
    :type step: dict
    :return: Tuple of (activation CV term name, electron energy in eV or None, reaction time in ms or None). The
        activation is ExD if an electron energy is set, ETD if only a reaction time is set, and CID otherwise.
    :rtype: tuple
    """
    electron_energy = step.get('ElectronEnergy')
    reaction_time = step.get('ReactionTime')
    if electron_energy is not None and not pd.isna(electron_energy):
        return 'electron activated dissociation', float(electron_energy), \
            None if reaction_time is None or pd.isna(reaction_time) else float(reaction_time)
    if reaction_time is not None and not pd.isna(reaction_time):
        return 'electron transfer dissociation', None, float(reaction_time)
    return 'collision-induced dissociation', None, None
