"""Maps resolved trial variables (t-vars) to GratingInstance parameters."""

from grating_renderer import GratingInstance


class StimulusManager:
    """Translates trial variables from ExperimentModel into GratingInstances.

    Owns the 4 grating instances and updates them from the model's t-variables
    at the start of each trial (or when the protocol engine signals a state change).
    """

    def __init__(self):
        self.stim_one = GratingInstance("stimOne")
        self.mask_one = GratingInstance("maskOne")
        self.stim_two = GratingInstance("stimTwo")
        self.mask_two = GratingInstance("maskTwo")
        self.all_gratings = [self.stim_one, self.mask_one, self.stim_two, self.mask_two]

    def apply_trial_params(self, model):
        """Copy t-variables from model into grating instances.

        Called once at the start of each trial, after the protocol engine
        has resolved all randomization into t-variables.
        """
        m = model

        # Stimulus One
        self.stim_one.direction_deg = m.get_float('tStimOneGratingDirectionDeg')
        self.stim_one.contrast = m.get_float('tStimOneGratingContrast')
        self.stim_one.spatial_freq = m.get_float('tStimOneGratingSpatialFreqCPD', 0.05)
        self.stim_one.temporal_freq = m.get_float('tStimOneGratingTemporalFreqCPS')
        self.stim_one.phase_deg = m.get_float('tStimOneGratingPhaseDeg')
        self.stim_one.diameter_deg = m.get_float('tStimOneGratingDiameterDeg', 200.0)
        self.stim_one.azimuth_deg = m.get_float('tStimOneGratingAzimuthDeg')
        self.stim_one.elevation_deg = m.get_float('tStimOneGratingElevationDeg')
        self.stim_one.parent_spatial_freq = None  # stim uses own SF

        # Mask One: uses parent stim's SF for spatial_frequency and speed
        self.mask_one.direction_deg = m.get_float('tMaskOneGratingDirectionDeg')
        self.mask_one.contrast = m.get_float('tMaskOneGratingContrast')
        self.mask_one.phase_deg = m.get_float('tMaskOneGratingPhaseDeg')
        self.mask_one.temporal_freq = m.get_float('tMaskOneGratingTemporalFreqCPS')
        # Mask shares parent stim's SF, diameter, position
        self.mask_one.spatial_freq = m.get_float('tMaskOneGratingTemporalFreqCPS')  # own TF for speed calc
        self.mask_one.parent_spatial_freq = m.get_float('tStimOneGratingSpatialFreqCPD', 0.05)
        self.mask_one.diameter_deg = m.get_float('tStimOneGratingDiameterDeg', 200.0)
        self.mask_one.azimuth_deg = m.get_float('tStimOneGratingAzimuthDeg')
        self.mask_one.elevation_deg = m.get_float('tStimOneGratingElevationDeg')

        # Stimulus Two
        self.stim_two.direction_deg = m.get_float('tStimTwoGratingDirectionDeg')
        self.stim_two.contrast = m.get_float('tStimTwoGratingContrast')
        self.stim_two.spatial_freq = m.get_float('tStimTwoGratingSpatialFreqCPD', 0.05)
        self.stim_two.temporal_freq = m.get_float('tStimTwoGratingTemporalFreqCPS')
        self.stim_two.phase_deg = m.get_float('tStimTwoGratingPhaseDeg')
        self.stim_two.diameter_deg = m.get_float('tStimTwoGratingDiameterDeg', 200.0)
        self.stim_two.azimuth_deg = m.get_float('tStimTwoGratingAzimuthDeg')
        self.stim_two.elevation_deg = m.get_float('tStimTwoGratingElevationDeg')
        self.stim_two.parent_spatial_freq = None

        # Mask Two: uses parent stim Two's SF
        self.mask_two.direction_deg = m.get_float('tMaskTwoGratingDirectionDeg')
        self.mask_two.contrast = m.get_float('tMaskTwoGratingContrast')
        self.mask_two.phase_deg = m.get_float('tMaskTwoGratingPhaseDeg')
        self.mask_two.temporal_freq = m.get_float('tMaskTwoGratingTemporalFreqCPS')
        self.mask_two.spatial_freq = m.get_float('tMaskTwoGratingTemporalFreqCPS')
        self.mask_two.parent_spatial_freq = m.get_float('tStimTwoGratingSpatialFreqCPD', 0.05)
        self.mask_two.diameter_deg = m.get_float('tStimTwoGratingDiameterDeg', 200.0)
        self.mask_two.azimuth_deg = m.get_float('tStimTwoGratingAzimuthDeg')
        self.mask_two.elevation_deg = m.get_float('tStimTwoGratingElevationDeg')

    def apply_gaussian_params(self, model):
        """Apply shared gaussian mask params (gratingStd, gratingMean, gratingEdge)."""
        std = model.get_float('gratingStd', 0.3)
        mean = model.get_float('gratingMean', 0.1)
        edge = model.get_float('gratingEdge', 0.125)
        for g in self.all_gratings:
            g.std_dev = std
            g.mean = mean
            g.edge_width = edge

    def get_active_gratings_for_state(self, state: str, model) -> list[GratingInstance]:
        """Return which gratings should be visible for the current protocol state.

        States: STIM_ONE_ON, STIM_TWO_ON, etc.
        """
        do_two_together = model.get_bool('doTwoStimTogether')

        if state == "STIM_ONE_ON" or state == "ROTATE":
            active = []
            if self.stim_one.contrast > 0:
                active.append(self.stim_one)
            if self.mask_one.contrast > 0:
                active.append(self.mask_one)
            if do_two_together:
                if self.stim_two.contrast > 0:
                    active.append(self.stim_two)
                if self.mask_two.contrast > 0:
                    active.append(self.mask_two)
            return active

        elif state == "STIM_TWO_ON":
            if do_two_together:
                return []  # already shown with stim one
            active = []
            if self.stim_two.contrast > 0:
                active.append(self.stim_two)
            if self.mask_two.contrast > 0:
                active.append(self.mask_two)
            return active

        else:
            # ITI, ISI, END_TRIAL, IDLE — no gratings
            return []

    def reset_phases(self):
        """Reset phase accumulators on all gratings."""
        for g in self.all_gratings:
            g.reset_phase()

    def update_phases(self, dt: float):
        """Advance phase on all gratings."""
        for g in self.all_gratings:
            g.update_phase(dt)
