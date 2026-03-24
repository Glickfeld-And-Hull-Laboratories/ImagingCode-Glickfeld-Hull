"""MWorks protocol state machine — drives trial-by-trial execution.

Translates the JuiceOnHoldProtocol from twoStim_2P_Frames.mwel into Python.

State machine:
  IDLE → INITIALIZATION → ITI_WAIT → STIM_ONE_ON → [ROTATE] → STIM_ONE_OFF →
         STIM_TWO_ON → STIM_TWO_OFF → END_TRIAL → INITIALIZATION → ...
"""

import math
import random

from experiment_model import ExperimentModel
from selection_variable import SelectionVariable
from stimulus_manager import StimulusManager
from grating_renderer import GratingInstance


class ProtocolEngine:
    """Runs the MWorks trial state machine at frameRateHz."""

    # State names
    IDLE = "IDLE"
    INITIALIZATION = "INITIALIZATION"
    ITI_WAIT = "ITI_WAIT"
    STIM_ONE_ON = "STIM_ONE_ON"
    ROTATE = "ROTATE"
    STIM_ONE_OFF = "STIM_ONE_OFF"
    STIM_TWO_ON = "STIM_TWO_ON"
    STIM_TWO_OFF = "STIM_TWO_OFF"
    END_TRIAL = "END_TRIAL"
    PAUSED = "PAUSED"

    def __init__(self, model: ExperimentModel, stim_mgr: StimulusManager):
        self.model = model
        self.stim_mgr = stim_mgr

        # Selection variables (4 independent pools of 80)
        self.sv_stim = SelectionVariable("svStimNumber")
        self.sv_stim2 = SelectionVariable("svStimNumber2")
        self.sv_mask = SelectionVariable("svMaskNumber")
        self.sv_mask2 = SelectionVariable("svMaskNumber2")

        # Counters matching MWEL
        self.n_stim_accepted = 0
        self.n_stim_accepted2 = 0
        self.n_mask_accepted = 0
        self.n_mask_accepted2 = 0

        # State
        self.state = self.IDLE
        self._pre_pause_state = self.IDLE
        self.counter = 0           # frame counter (matches MWEL 'counter')
        self.ntrials = 0
        self.trial_number = 0      # total trials completed
        self.is_first_trial = True

        # Per-trial timing anchors (frame numbers)
        self._c_iti_start = 0
        self._c_stim_one_on = 0
        self._c_stim_one_off = 0
        self._c_stim_two_on = 0
        self._c_stim_two_off = 0
        self._c_now = 0

        # Per-trial computed frame counts
        self._n_iti_wait_frames = 0
        self._n_stim_one_frames_on = 0
        self._n_stim_two_frames_on = 0
        self._n_frames_isi = 0
        self._n_deg_rotation_per_frame = 1.0

        # Trial info for GUI
        self.last_trial_time_ms = 0
        self._this_trial_start_ms = -1
        self._last_trial_start_ms = -1

    # ── Public API ──

    def start(self):
        """Start the protocol from the beginning."""
        self._reset_selections()
        self.counter = 0
        self.ntrials = 0
        self.trial_number = 0
        self.is_first_trial = True
        self.state = self.INITIALIZATION

    def stop(self):
        """Stop the protocol."""
        self.state = self.IDLE

    def pause(self):
        """Pause the protocol."""
        if self.state != self.IDLE and self.state != self.PAUSED:
            self._pre_pause_state = self.state
            self.state = self.PAUSED

    def resume(self):
        """Resume from pause."""
        if self.state == self.PAUSED:
            self.state = self._pre_pause_state

    def skip_trial(self):
        """Skip to the next trial's initialization."""
        if self.state not in (self.IDLE, self.PAUSED):
            self._do_end_trial()
            self.state = self.INITIALIZATION

    def tick(self) -> list[GratingInstance]:
        """Advance one frame. Returns list of active gratings to render.

        Called at frameRateHz from main loop's frame accumulator.
        """
        if self.state == self.IDLE or self.state == self.PAUSED:
            return []

        self.counter += 1
        active: list[GratingInstance] = []

        if self.state == self.INITIALIZATION:
            self._do_initialization()
            self.state = self.ITI_WAIT
            self._c_iti_start = self.counter

        elif self.state == self.ITI_WAIT:
            # Wait for ITI frames
            if self.counter >= self._c_iti_start + self._n_iti_wait_frames:
                if self.is_first_trial and self.counter <= 19:
                    pass  # MWEL: wait for counter > 19 on first trial
                else:
                    self.state = self.STIM_ONE_ON
                    self._enter_stim_one_on()

        elif self.state == self.STIM_ONE_ON:
            active = self.stim_mgr.get_active_gratings_for_state(self.STIM_ONE_ON, self.model)
            # Check for rotation
            do_rotate = self.model.get_bool('tDoRotate')
            if do_rotate and self.counter >= self._c_stim_one_on + 1:
                self.state = self.ROTATE
                self._c_now = self.counter
            elif self.counter >= self._c_stim_one_on + self._n_stim_one_frames_on:
                self.state = self.STIM_ONE_OFF
                self._c_stim_one_off = self.counter

        elif self.state == self.ROTATE:
            # Update direction each frame
            frames_since = self.counter - self._c_stim_one_on
            base_dir = self.model.get_float('tStimOneGratingDirectionDeg')
            new_dir = base_dir + frames_since * self._n_deg_rotation_per_frame
            self.model.set('setStimOneGratingDirectionDeg', new_dir)
            self.stim_mgr.stim_one.direction_deg = new_dir

            if self.model.get_bool('tDoMaskRotate'):
                mask_offset = self.model.get_float('maskOneGratingDirectionDeg') - self.model.get_float('stimOneGratingDirectionDeg')
                self.stim_mgr.mask_one.direction_deg = new_dir + mask_offset

            active = self.stim_mgr.get_active_gratings_for_state(self.ROTATE, self.model)

            if self.counter >= self._c_stim_one_on + self._n_stim_one_frames_on:
                self.state = self.STIM_ONE_OFF
                self._c_stim_one_off = self.counter
            else:
                self._c_now = self.counter

        elif self.state == self.STIM_ONE_OFF:
            # ISI period — no gratings
            if self.counter >= self._c_stim_one_off + self._n_frames_isi:
                self.state = self.STIM_TWO_ON
                self._enter_stim_two_on()

        elif self.state == self.STIM_TWO_ON:
            active = self.stim_mgr.get_active_gratings_for_state(self.STIM_TWO_ON, self.model)
            if self.counter >= self._c_stim_two_on + self._n_stim_two_frames_on:
                self.state = self.STIM_TWO_OFF
                self._c_stim_two_off = self.counter

        elif self.state == self.STIM_TWO_OFF:
            # Wait 3 frames then end trial (matches MWEL)
            if self.counter >= self._c_stim_two_off + 3:
                self._do_end_trial()
                self.state = self.INITIALIZATION

        return active

    def get_trial_info(self) -> dict:
        """Return current trial info for GUI display."""
        m = self.model
        return {
            'state': self.state,
            'trial_number': self.trial_number,
            'counter': self.counter,
            'stim_number': m.get_int('tStimulusNumber'),
            'stim_number2': m.get_int('tStimulusNumber2'),
            'mask_number': m.get_int('tMaskNumber'),
            'block2': m.get_bool('tBlock2TrialNumber'),
            'stim1_dir': m.get_float('tStimOneGratingDirectionDeg'),
            'stim1_con': m.get_float('tStimOneGratingContrast'),
            'mask1_dir': m.get_float('tMaskOneGratingDirectionDeg'),
            'mask1_con': m.get_float('tMaskOneGratingContrast'),
            'stim2_dir': m.get_float('tStimTwoGratingDirectionDeg'),
            'stim2_con': m.get_float('tStimTwoGratingContrast'),
            'do_mask': m.get_bool('tDoMask'),
            'isi_ms': m.get_float('tISITimeMs'),
            'stim1_on_ms': m.get_float('tStimOneGratingOnTimeMs'),
            'n_stim_accepted': self.n_stim_accepted,
        }

    # ── Private: state entry actions ──

    def _enter_stim_one_on(self):
        """Actions on entering StimOneOn state."""
        self.ntrials += 1
        self.trial_number += 1
        self._c_stim_one_on = self.counter
        # Reset grating phases for new stimulus presentation
        self.stim_mgr.reset_phases()

    def _enter_stim_two_on(self):
        """Actions on entering StimTwoOn state."""
        self._c_stim_two_on = self.counter
        # Reset stim two phases
        self.stim_mgr.stim_two.reset_phase()
        self.stim_mgr.mask_two.reset_phase()

    # ── Private: initialization (all randomization happens here) ──

    def _do_initialization(self):
        """Resolve all trial variables — matches MWEL Initialization state."""
        m = self.model

        if self.ntrials > 0:
            self.is_first_trial = False

        # Draw from selection variables
        stim_num = self.sv_stim.current_value
        stim_num2 = self.sv_stim2.current_value
        mask_num = self.sv_mask.current_value
        mask_num2 = self.sv_mask2.current_value

        m.set('tStimulusNumber', stim_num)
        m.set('tStimulusNumber2', stim_num2)
        m.set('tMaskNumber', mask_num)
        m.set('tMaskNumber2', mask_num2)

        # Block 2 determination
        do_block2 = m.get_bool('doBlock2')
        if do_block2:
            if not m.get_bool('doBlock2SeparateOdds'):
                m.set('tBlock2TrialNumber', 1 if stim_num <= 40 else 0)
            else:
                m.set('tBlock2TrialNumber', 1 if stim_num <= m.get_float('block2TrPer80') else 0)
        else:
            m.set('tBlock2TrialNumber', 0)

        is_b2 = m.get_bool('tBlock2TrialNumber')

        # Set default t-variables from base parameters
        self._set_default_tvars(m)

        # Handle matchStimOneParameters
        if m.get_bool('matchStimOneParameters'):
            self._apply_match_stim_one(m)
        else:
            self._apply_stim_two_defaults(m)

        # Block 2 overrides
        if do_block2 and is_b2:
            self._apply_block2_overrides(m)

        # Single condition randomization (not matrix)
        if not m.get_bool('tDoMatrix'):
            self._apply_single_randomization(m, stim_num, stim_num2, mask_num, mask_num2, is_b2)
        else:
            self._apply_matrix_randomization(m, stim_num, stim_num2, is_b2)

        # Post-randomization fixups
        self._post_randomization_fixups(m, stim_num, stim_num2, is_b2)

        # Frame conversion
        frame_rate = m.get_float('frameRateHz', 30.0)
        stim1_on_ms = m.get_float('tStimOneGratingOnTimeMs', 500.0)
        stim2_on_ms = m.get_float('tStimTwoGratingOnTimeMs', 500.0)
        isi_ms = m.get_float('tISITimeMs', 250.0)
        rotate_speed = m.get_float('gratingRotateSpeedDPS', 20.0)

        self._n_stim_one_frames_on = math.ceil((stim1_on_ms / 1000.0) * frame_rate)
        self._n_stim_two_frames_on = math.ceil((stim2_on_ms / 1000.0) * frame_rate)
        self._n_frames_isi = math.ceil((isi_ms / 1000.0) * frame_rate)
        self._n_deg_rotation_per_frame = rotate_speed / frame_rate

        if m.get_bool('doRotate'):
            self._n_stim_one_frames_on = math.ceil(360.0 / self._n_deg_rotation_per_frame)

        # ITI calculation
        iti_ms = m.get_float('itiTimeMs', 5000.0)
        m.set('tItiWaitTimeMs', iti_ms)
        self._n_iti_wait_frames = math.ceil((iti_ms / 1000.0) * frame_rate)

        # Apply trial params to grating instances
        self.stim_mgr.apply_trial_params(m)
        self.stim_mgr.apply_gaussian_params(m)

    def _set_default_tvars(self, m):
        """Copy base parameters to t-variables."""
        m.set('tISITimeMs', m.get_float('isiTimeMs', 250.0))
        m.set('tStimOneDoVisualStim', m.get_float('stimOneDoVisualStim', 1.0))
        m.set('tStimOneDoAuditoryStim', m.get_float('stimOneDoAuditoryStim', 0.0))
        m.set('tStimOneGratingOnTimeMs', m.get_float('stimOneGratingOnTimeMs', 500.0))
        m.set('tStimOneGratingDirectionDeg', m.get_float('stimOneGratingDirectionDeg'))
        m.set('tStimOneGratingContrast', m.get_float('stimOneGratingContrast', 1.0))
        m.set('tStimOneGratingElevationDeg', m.get_float('stimOneGratingElevationDeg'))
        m.set('tStimOneGratingAzimuthDeg', m.get_float('stimOneGratingAzimuthDeg'))
        m.set('tStimOneGratingDiameterDeg', m.get_float('stimOneGratingDiameterDeg', 200.0))
        m.set('tStimOneGratingSpatialFreqCPD', m.get_float('stimOneGratingSpatialFreqCPD', 0.05))
        m.set('tStimOneGratingTemporalFreqCPS', m.get_float('stimOneGratingTemporalFreqCPS'))
        m.set('tStimOneGratingPhaseDeg', m.get_float('stimOneGratingPhaseDeg'))
        m.set('tMaskOneGratingDirectionDeg', m.get_float('maskOneGratingDirectionDeg'))
        m.set('tMaskOneGratingContrast', m.get_float('maskOneGratingContrast'))
        m.set('tMaskOneGratingPhaseDeg', m.get_float('maskOneGratingPhaseDeg'))
        m.set('tMaskOneGratingTemporalFreqCPS', m.get_float('maskOneGratingTemporalFreqCPS'))
        m.set('tStimTwoDoVisualStim', m.get_float('stimTwoDoVisualStim', 1.0))
        m.set('tStimTwoDoAuditoryStim', m.get_float('stimTwoDoAuditoryStim'))
        m.set('tDoMask', m.get_bool('doMask'))
        m.set('tDoStimTwoMask', m.get_bool('doStimTwoMask'))
        m.set('tDoMatrix', m.get_bool('doMatrix'))
        m.set('tDoRandISITime', m.get_bool('doRandISITime'))
        m.set('tDoRandStimOnTime', m.get_bool('doRandStimOnTime'))
        m.set('tDoRandCon', m.get_bool('doRandCon'))
        m.set('tDoRandDir', m.get_bool('doRandDir'))
        m.set('tDoStimTwoRandDir', m.get_bool('doStimTwoRandDir'))
        m.set('tDoRandMaskCon', m.get_bool('doRandMaskCon'))
        m.set('tDoRandMaskDir', m.get_bool('doRandMaskDir'))
        m.set('tDoRandMaskPhase', m.get_bool('doRandMaskPhase'))
        m.set('tDoRandMaskTF', m.get_bool('doRandMaskTF'))
        m.set('tDoRandPos', m.get_bool('doRandPos'))
        m.set('tDoRandDiameter', m.get_bool('doRandDiameter'))
        m.set('tDoRandSF', m.get_bool('doRandSF'))
        m.set('tDoRandTF', m.get_bool('doRandTF'))
        m.set('tDoRandPhase', m.get_bool('doRandPhase'))
        m.set('tDoRotate', m.get_bool('doRotate'))
        m.set('tDoMaskRotate', m.get_bool('doMaskRotate'))
        m.set('tDoRandSoundAmp', m.get_bool('doRandSoundAmp'))
        m.set('tDoRandSoundFreq', m.get_bool('doRandSoundFreq'))
        m.set('tStimOneSoundAmplitude', m.get_float('stimOneSoundAmplitude'))

    def _apply_match_stim_one(self, m):
        """matchStimOneParameters: copy stim one params to stim two."""
        m.set('tStimTwoGratingOnTimeMs', m.get_float('stimOneGratingOnTimeMs', 500.0))
        m.set('tStimTwoGratingDirectionDeg', m.get_float('stimOneGratingDirectionDeg'))
        m.set('tStimTwoGratingContrast', m.get_float('stimOneGratingContrast', 1.0))
        m.set('tStimTwoGratingElevationDeg', m.get_float('stimOneGratingElevationDeg'))
        m.set('tStimTwoGratingAzimuthDeg', m.get_float('stimOneGratingAzimuthDeg'))
        m.set('tStimTwoGratingDiameterDeg', m.get_float('stimOneGratingDiameterDeg', 200.0))
        m.set('tStimTwoGratingSpatialFreqCPD', m.get_float('stimOneGratingSpatialFreqCPD', 0.05))
        m.set('tStimTwoGratingTemporalFreqCPS', m.get_float('stimOneGratingTemporalFreqCPS'))
        m.set('tStimTwoGratingPhaseDeg', m.get_float('stimOneGratingPhaseDeg'))
        m.set('tMaskTwoGratingDirectionDeg', m.get_float('maskOneGratingDirectionDeg'))
        m.set('tMaskTwoGratingContrast', m.get_float('maskOneGratingContrast'))
        m.set('tMaskTwoGratingPhaseDeg', m.get_float('maskOneGratingPhaseDeg'))
        m.set('tStimTwoSoundAmplitude', m.get_float('stimOneSoundAmplitude'))

    def _apply_stim_two_defaults(self, m):
        """Set stim two t-vars from stim two base params."""
        m.set('tStimTwoGratingOnTimeMs', m.get_float('stimTwoGratingOnTimeMs', 500.0))
        m.set('tStimTwoGratingDirectionDeg', m.get_float('stimTwoGratingDirectionDeg'))
        m.set('tStimTwoGratingContrast', m.get_float('stimTwoGratingContrast', 1.0))
        m.set('tStimTwoGratingElevationDeg', m.get_float('stimTwoGratingElevationDeg'))
        m.set('tStimTwoGratingAzimuthDeg', m.get_float('stimTwoGratingAzimuthDeg'))
        m.set('tStimTwoGratingDiameterDeg', m.get_float('stimTwoGratingDiameterDeg', 200.0))
        m.set('tStimTwoGratingSpatialFreqCPD', m.get_float('stimTwoGratingSpatialFreqCPD', 0.05))
        m.set('tStimTwoGratingTemporalFreqCPS', m.get_float('stimTwoGratingTemporalFreqCPS'))
        m.set('tStimTwoGratingPhaseDeg', m.get_float('stimTwoGratingPhaseDeg'))
        m.set('tMaskTwoGratingDirectionDeg', m.get_float('maskTwoGratingDirectionDeg'))
        m.set('tMaskTwoGratingContrast', m.get_float('maskTwoGratingContrast'))
        m.set('tMaskTwoGratingPhaseDeg', m.get_float('maskTwoGratingPhaseDeg'))
        m.set('tStimTwoSoundAmplitude', m.get_float('stimTwoSoundAmplitude'))
        m.set('tMaskTwoGratingTemporalFreqCPS', m.get_float('maskTwoGratingTemporalFreqCPS'))

    def _apply_block2_overrides(self, m):
        """Apply Block 2 parameter overrides when tBlock2TrialNumber is true."""
        if not m.get_bool('block2MatchB1VisStim'):
            m.set('tDoMask', m.get_bool('block2DoMask'))
            m.set('tDoStimTwoMask', m.get_bool('block2DoStimTwoMask'))
            m.set('tDoMatrix', m.get_bool('block2DoMatrix'))
            m.set('tDoRandISITime', m.get_bool('block2DoRandISITime'))
            m.set('tDoRandStimOnTime', m.get_bool('block2DoRandStimOnTime'))
            m.set('tDoRandCon', m.get_bool('block2DoRandCon'))
            m.set('tDoRandDir', m.get_bool('block2DoRandDir'))
            m.set('tDoStimTwoRandDir', m.get_bool('block2DoStimTwoRandDir'))
            m.set('tDoRandPos', m.get_bool('block2DoRandPos'))
            m.set('tDoRandDiameter', m.get_bool('block2DoRandDiameter'))
            m.set('tDoRandSF', m.get_bool('block2DoRandSF'))
            m.set('tDoRandTF', m.get_bool('block2DoRandTF'))
            m.set('tDoRandPhase', m.get_bool('block2DoRandPhase'))
            m.set('tDoRotate', m.get_bool('block2DoRotate'))
            m.set('tDoMaskRotate', m.get_bool('block2DoMaskRotate'))
            m.set('tDoRandMaskCon', m.get_bool('block2DoRandMaskCon'))
            m.set('tDoRandMaskDir', m.get_bool('block2DoRandMaskDir'))
            m.set('tDoRandMaskPhase', m.get_bool('block2DoRandMaskPhase'))
            m.set('tDoRandSoundAmp', m.get_bool('block2DoRandSoundAmp'))
            m.set('tDoRandSoundFreq', m.get_bool('block2DoRandSoundFreq'))
            m.set('tISITimeMs', m.get_float('block2ISITimeMs', 250.0))
            m.set('tStimOneDoVisualStim', m.get_float('block2StimOneDoVisualStim'))
            m.set('tStimTwoDoVisualStim', m.get_float('block2StimTwoDoVisualStim'))
            m.set('tStimOneGratingOnTimeMs', m.get_float('block2StimOneGratingOnTimeMs'))
            m.set('tStimOneGratingDirectionDeg', m.get_float('block2StimOneGratingDirectionDeg'))
            m.set('tStimOneGratingContrast', m.get_float('block2StimOneGratingContrast'))
            m.set('tStimOneGratingElevationDeg', m.get_float('block2StimOneGratingElevationDeg'))
            m.set('tStimOneGratingAzimuthDeg', m.get_float('block2StimOneGratingAzimuthDeg'))
            m.set('tStimOneGratingDiameterDeg', m.get_float('block2StimOneGratingDiameterDeg'))
            m.set('tStimOneGratingSpatialFreqCPD', m.get_float('block2StimOneGratingSpatialFreqCPD'))
            m.set('tStimOneGratingTemporalFreqCPS', m.get_float('block2StimOneGratingTemporalFreqCPS'))
            m.set('tStimOneGratingPhaseDeg', m.get_float('block2StimOneGratingPhaseDeg'))
            m.set('tMaskOneGratingDirectionDeg', m.get_float('block2MaskOneGratingDirectionDeg'))
            m.set('tMaskOneGratingContrast', m.get_float('block2MaskOneGratingContrast'))
            m.set('tMaskOneGratingPhaseDeg', m.get_float('block2MaskOneGratingPhaseDeg'))

            if m.get_bool('block2MatchStimOneParameters'):
                m.set('tStimTwoGratingDirectionDeg', m.get_float('block2StimOneGratingDirectionDeg'))
                m.set('tStimTwoGratingContrast', m.get_float('block2StimOneGratingContrast'))
                m.set('tStimTwoGratingElevationDeg', m.get_float('block2StimOneGratingElevationDeg'))
                m.set('tStimTwoGratingAzimuthDeg', m.get_float('block2StimOneGratingAzimuthDeg'))
                m.set('tStimTwoGratingDiameterDeg', m.get_float('block2StimOneGratingDiameterDeg'))
                m.set('tStimTwoGratingSpatialFreqCPD', m.get_float('block2StimOneGratingSpatialFreqCPD'))
                m.set('tStimTwoGratingTemporalFreqCPS', m.get_float('block2StimOneGratingTemporalFreqCPS'))
                m.set('tStimTwoGratingPhaseDeg', m.get_float('block2StimOneGratingPhaseDeg'))
                m.set('tMaskTwoGratingDirectionDeg', m.get_float('block2MaskOneGratingDirectionDeg'))
                m.set('tMaskTwoGratingContrast', m.get_float('block2MaskOneGratingContrast'))
                m.set('tMaskTwoGratingPhaseDeg', m.get_float('block2MaskOneGratingPhaseDeg'))
            else:
                m.set('tStimTwoGratingOnTimeMs', m.get_float('block2StimTwoGratingOnTimeMs'))
                m.set('tStimTwoGratingDirectionDeg', m.get_float('block2StimTwoGratingDirectionDeg'))
                m.set('tStimTwoGratingContrast', m.get_float('block2StimTwoGratingContrast'))
                m.set('tStimTwoGratingElevationDeg', m.get_float('block2StimTwoGratingElevationDeg'))
                m.set('tStimTwoGratingAzimuthDeg', m.get_float('block2StimTwoGratingAzimuthDeg'))
                m.set('tStimTwoGratingDiameterDeg', m.get_float('block2StimTwoGratingDiameterDeg'))
                m.set('tStimTwoGratingSpatialFreqCPD', m.get_float('block2StimTwoGratingSpatialFreqCPD'))
                m.set('tStimTwoGratingTemporalFreqCPS', m.get_float('block2StimTwoGratingTemporalFreqCPS'))
                m.set('tStimTwoGratingPhaseDeg', m.get_float('block2StimTwoGratingPhaseDeg'))
                m.set('tMaskTwoGratingDirectionDeg', m.get_float('block2MaskTwoGratingDirectionDeg'))
                m.set('tMaskTwoGratingContrast', m.get_float('block2MaskTwoGratingContrast'))
                m.set('tMaskTwoGratingPhaseDeg', m.get_float('block2MaskTwoGratingPhaseDeg'))

    # ── Randomization helpers ──

    @staticmethod
    def _linear_rand(base, step, n_steps, index):
        """Linear randomization: base + step * (index % n_steps)."""
        if n_steps <= 0:
            return base
        return base + step * (index % n_steps)

    @staticmethod
    def _log_rand(base, step_log, n_steps, index):
        """Logarithmic randomization: base * pow(step_log, index % n_steps)."""
        if n_steps <= 0:
            return base
        return base * pow(step_log, index % n_steps)

    def _apply_single_randomization(self, m, sn, sn2, mn, mn2, is_b2):
        """Apply single-condition randomization (non-matrix mode)."""
        match_s1 = m.get_bool('matchStimOneParameters')
        do_b2 = m.get_bool('doBlock2')
        b2_match_vis = m.get_bool('block2MatchB1VisStim')
        b2_match_s1 = m.get_bool('block2MatchStimOneParameters')

        # ISI time randomization
        if m.get_bool('tDoRandISITime'):
            t = sn % max(1, int(m.get_float('isiTimeN', 5)))
            m.set('tISITimeMs', m.get_float('isiTimeMs', 250) * pow(m.get_float('isiTimeStepLog', 2), t))
            if do_b2 and is_b2 and not b2_match_vis:
                t = sn % max(1, int(m.get_float('block2ISITimeN', 5)))
                m.set('tISITimeMs', m.get_float('block2ISITimeMs', 250) + pow(m.get_float('block2ISITimeStepLog', 2), t))

        # Stim on time randomization
        if m.get_bool('tDoRandStimOnTime') and not is_b2:
            t1 = sn % max(1, int(m.get_float('stimOneGratingOnTimeStepN', 5)))
            m.set('tStimOneGratingOnTimeMs', m.get_float('stimOneGratingOnTimeMs', 500) * pow(m.get_float('stimOneGratingOnTimeStepLog', 2), t1))
            if match_s1:
                m.set('tStimTwoGratingOnTimeMs', m.get_float('tStimOneGratingOnTimeMs'))
            else:
                t2 = sn2 % max(1, int(m.get_float('stimTwoGratingOnTimeStepN', 5)))
                m.set('tStimTwoGratingOnTimeMs', m.get_float('stimTwoGratingOnTimeMs', 500) * pow(m.get_float('stimTwoGratingOnTimeStepLog', 2), t2))

        # Contrast randomization
        if m.get_bool('tDoRandCon'):
            t1 = sn % max(1, int(m.get_float('stimOneGratingContrastStepN', 5)))
            c1 = m.get_float('stimOneGratingContrast', 1) * pow(m.get_float('stimOneGratingContrastStepLog', 2), t1)
            if m.get_bool('doZeroCon') and sn < m.get_int('zeroConPer80'):
                c1 = 0
            m.set('tStimOneGratingContrast', c1)
            if match_s1:
                m.set('tStimTwoGratingContrast', c1)
            else:
                t2 = sn2 % max(1, int(m.get_float('stimTwoGratingContrastStepN', 5)))
                c2 = m.get_float('stimTwoGratingContrast', 1) * pow(m.get_float('stimTwoGratingContrastStepLog', 2), t2)
                if m.get_bool('doZeroCon') and sn2 < m.get_int('zeroConPer80'):
                    c2 = 0
                m.set('tStimTwoGratingContrast', c2)
            # Block 2 contrast override
            if do_b2 and is_b2 and not b2_match_vis:
                t1 = sn % max(1, int(m.get_float('block2StimOneGratingContrastStepN', 5)))
                c1 = m.get_float('block2StimOneGratingContrast', 1) * pow(m.get_float('block2StimOneGratingContrastStepLog', 2), t1)
                if m.get_bool('block2DoZeroCon') and sn < m.get_int('block2ZeroConPer80'):
                    c1 = 0
                m.set('tStimOneGratingContrast', c1)
                if b2_match_s1:
                    m.set('tStimTwoGratingContrast', c1)
                else:
                    t2 = sn2 % max(1, int(m.get_float('block2StimTwoGratingContrastStepN', 5)))
                    c2 = m.get_float('block2StimTwoGratingContrast', 1) * pow(m.get_float('block2StimTwoGratingContrastStepLog', 2), t2)
                    if m.get_bool('block2DoZeroCon') and sn2 < m.get_int('block2ZeroConPer80'):
                        c2 = 0
                    m.set('tStimTwoGratingContrast', c2)

        # Zero contrast catch trials (no doRandCon case)
        if m.get_bool('doZeroCon') and not m.get_bool('doRandCon'):
            if (is_b2 and sn < m.get_int('zeroConPer80')) or (not is_b2 and sn > 80 - m.get_int('zeroConPer80')):
                m.set('tStimOneGratingContrast', 0)

        # Direction randomization
        if m.get_bool('tDoRandDir'):
            n = max(1, int(m.get_float('stimOneGratingDirectionStepN', 12)))
            t1 = sn % n
            d1 = m.get_float('stimOneGratingDirectionDeg') + m.get_float('stimOneGratingDirectionStepDeg', 30) * t1
            m.set('tStimOneGratingDirectionDeg', d1)
            m.set('tMaskOneGratingDirectionDeg', d1 + m.get_float('maskOneGratingDirectionDeg'))
            if match_s1:
                m.set('tStimTwoGratingDirectionDeg', d1)
            else:
                n2 = max(1, int(m.get_float('stimTwoGratingDirectionStepN', 12)))
                t2 = sn2 % n2
                d2 = m.get_float('stimTwoGratingDirectionDeg') + m.get_float('stimTwoGratingDirectionStepDeg', 30) * t2
                m.set('tStimTwoGratingDirectionDeg', d2)
            # Block 2
            if do_b2 and is_b2 and not b2_match_vis:
                n = max(1, int(m.get_float('block2StimOneGratingDirectionStepN', 12)))
                t1 = sn % n
                d1 = m.get_float('block2StimOneGratingDirectionDeg') + m.get_float('block2StimOneGratingDirectionStepDeg', 30) * t1
                m.set('tStimOneGratingDirectionDeg', d1)
                if b2_match_s1:
                    m.set('tStimTwoGratingDirectionDeg', d1)
                else:
                    n2 = max(1, int(m.get_float('block2StimTwoGratingDirectionStepN', 12)))
                    t2 = sn2 % n2
                    d2 = m.get_float('block2StimTwoGratingDirectionDeg') + m.get_float('block2StimTwoGratingDirectionStepDeg', 30) * t2
                    m.set('tStimTwoGratingDirectionDeg', d2)

        # Stim Two independent direction randomization
        if m.get_bool('tDoStimTwoRandDir'):
            n2 = max(1, int(m.get_float('stimTwoGratingDirectionStepN', 12)))
            t2 = sn2 % n2
            d2 = m.get_float('stimTwoGratingDirectionDeg') + m.get_float('stimTwoGratingDirectionStepDeg', 30) * t2
            m.set('tStimTwoGratingDirectionDeg', d2)
            m.set('tMaskTwoGratingDirectionDeg', d2 + m.get_float('maskTwoGratingDirectionDeg'))

        # Phase randomization
        if m.get_bool('tDoRandPhase'):
            t1 = sn % max(1, int(m.get_float('stimOneGratingPhaseStepN', 5)))
            m.set('tStimOneGratingPhaseDeg', m.get_float('stimOneGratingPhaseDeg') + m.get_float('stimOneGratingPhaseStepDeg', 2) * t1)
            if match_s1:
                m.set('tStimTwoGratingPhaseDeg', m.get_float('tStimOneGratingPhaseDeg'))
            else:
                t2 = sn2 % max(1, int(m.get_float('stimTwoGratingPhaseStepN', 5)))
                m.set('tStimTwoGratingPhaseDeg', m.get_float('stimTwoGratingPhaseDeg') + m.get_float('stimTwoGratingPhaseStepDeg', 2) * t2)
            if do_b2 and is_b2 and not b2_match_vis:
                t1 = sn % max(1, int(m.get_float('block2StimOneGratingPhaseStepN', 5)))
                m.set('tStimOneGratingPhaseDeg', m.get_float('block2StimOneGratingPhaseDeg') + m.get_float('block2StimOneGratingPhaseStepDeg', 2) * t1)
                if b2_match_s1:
                    m.set('tStimTwoGratingPhaseDeg', m.get_float('tStimOneGratingPhaseDeg'))
                else:
                    t2 = sn2 % max(1, int(m.get_float('block2StimTwoGratingPhaseStepN', 5)))
                    m.set('tStimTwoGratingPhaseDeg', m.get_float('block2StimTwoGratingPhaseDeg') + m.get_float('block2StimTwoGratingPhaseStepDeg', 2) * t2)

        # Diameter randomization
        if m.get_bool('tDoRandDiameter'):
            t1 = sn % max(1, int(m.get_float('stimOneGratingDiameterStepN', 5)))
            m.set('tStimOneGratingDiameterDeg', m.get_float('stimOneGratingDiameterDeg') * pow(m.get_float('stimOneGratingDiameterStepLog', 2), t1))
            if match_s1:
                m.set('tStimTwoGratingDiameterDeg', m.get_float('tStimOneGratingDiameterDeg'))
            else:
                t2 = sn2 % max(1, int(m.get_float('stimTwoGratingDiameterStepN', 5)))
                m.set('tStimTwoGratingDiameterDeg', m.get_float('stimTwoGratingDiameterDeg') * pow(m.get_float('stimTwoGratingDiameterStepLog', 2), t2))

        # TF randomization
        if m.get_bool('tDoRandTF'):
            t1 = sn % max(1, int(m.get_float('stimOneGratingTemporalFreqStepN', 5)))
            m.set('tStimOneGratingTemporalFreqCPS', m.get_float('stimOneGratingTemporalFreqCPS') * pow(m.get_float('stimOneGratingTemporalFreqStepLog', 2), t1))
            if match_s1:
                m.set('tStimTwoGratingTemporalFreqCPS', m.get_float('tStimOneGratingTemporalFreqCPS'))
            else:
                t2 = sn2 % max(1, int(m.get_float('stimTwoGratingTemporalFreqStepN', 5)))
                m.set('tStimTwoGratingTemporalFreqCPS', m.get_float('stimTwoGratingTemporalFreqCPS') * pow(m.get_float('stimTwoGratingTemporalFreqStepLog', 2), t2))

        # SF randomization
        if m.get_bool('tDoRandSF'):
            t1 = sn % max(1, int(m.get_float('stimOneGratingSpatialFreqStepN', 5)))
            m.set('tStimOneGratingSpatialFreqCPD', m.get_float('stimOneGratingSpatialFreqCPD') * pow(m.get_float('stimOneGratingSpatialFreqStepLog', 2), t1))
            if match_s1:
                m.set('tStimTwoGratingSpatialFreqCPD', m.get_float('tStimOneGratingSpatialFreqCPD'))
            else:
                t2 = sn2 % max(1, int(m.get_float('stimTwoGratingSpatialFreqStepN', 5)))
                m.set('tStimTwoGratingSpatialFreqCPD', m.get_float('stimTwoGratingSpatialFreqCPD') * pow(m.get_float('stimTwoGratingSpatialFreqStepLog', 2), t2))

        # Position randomization
        if m.get_bool('tDoRandPos'):
            self._randomize_position(m, sn, sn2, match_s1, is_b2, do_b2, b2_match_vis, b2_match_s1)

        # Mask randomization
        if m.get_bool('tDoMask'):
            self._randomize_mask(m, mn, mn2, sn, sn2, match_s1, is_b2, do_b2, b2_match_vis, b2_match_s1)

        # StimTwoMask
        if m.get_bool('tDoStimTwoMask'):
            mult_val = random.random()
            fract = m.get_float('fractMaskTrials')
            if mult_val >= fract:
                m.set('tDoStimTwoMask', 0)
                m.set('tMaskTwoGratingContrast', 0)
            else:
                m.set('tDoStimTwoMask', 1)
                m.set('tMaskTwoGratingContrast', m.get_float('stimTwoGratingContrast'))

    def _randomize_position(self, m, sn, sn2, match_s1, is_b2, do_b2, b2_match_vis, b2_match_s1):
        """Position (azimuth/elevation) randomization."""
        el_n = max(1, int(m.get_float('stimOneGratingElevationStepN', 12)))
        az_n = max(1, int(m.get_float('stimOneGratingAzimuthStepN', 12)))
        t1 = sn % (el_n * az_n)
        az_base = m.get_float('stimOneGratingAzimuthDeg')
        el_base = m.get_float('stimOneGratingElevationDeg')
        az_step = m.get_float('stimOneGratingAzimuthStepDeg', 30)
        el_step = m.get_float('stimOneGratingElevationStepDeg', 30)

        if t1 < az_n:
            m.set('tStimOneGratingAzimuthDeg', az_base + az_step * t1)
            m.set('tStimOneGratingElevationDeg', el_base)
        else:
            m.set('tStimOneGratingAzimuthDeg', az_base + az_step * (t1 % az_n))
            m.set('tStimOneGratingElevationDeg', el_base + el_step * (t1 // az_n))

        if match_s1:
            m.set('tStimTwoGratingAzimuthDeg', m.get_float('tStimOneGratingAzimuthDeg'))
            m.set('tStimTwoGratingElevationDeg', m.get_float('tStimOneGratingElevationDeg'))

    def _randomize_mask(self, m, mn, mn2, sn, sn2, match_s1, is_b2, do_b2, b2_match_vis, b2_match_s1):
        """Mask randomization: contrast, direction, phase, TF."""
        mult_val = random.random()
        fract = m.get_float('fractMaskTrials')
        if mult_val >= fract:
            m.set('tDoMask', 0)
            m.set('tMaskOneGratingContrast', 0)
            m.set('tMaskTwoGratingContrast', 0)
            return

        m.set('tDoMask', 1)
        do_rand_mask_con = m.get_bool('tDoRandMaskCon')
        do_rand_mask_dir = m.get_bool('tDoRandMaskDir')
        do_rand_mask_phase = m.get_bool('tDoRandMaskPhase')
        do_rand_mask_tf = m.get_bool('tDoRandMaskTF')

        # Mask contrast only (no dir)
        if do_rand_mask_con and not do_rand_mask_dir:
            n = max(1, int(m.get_float('maskOneGratingContrastStepN', 5)))
            t1 = mn % n
            c1 = m.get_float('maskOneGratingContrast') * pow(m.get_float('maskOneGratingContrastStepLog', 2), t1)
            if m.get_bool('doZeroCon') and mn < m.get_int('zeroConPer80'):
                c1 = 0
            m.set('tMaskOneGratingContrast', c1)
            if match_s1:
                m.set('tMaskTwoGratingContrast', c1)
            else:
                n2 = max(1, int(m.get_float('maskTwoGratingContrastStepN', 5)))
                t2 = mn2 % n2
                c2 = m.get_float('maskTwoGratingContrast') * pow(m.get_float('maskTwoGratingContrastStepLog', 2), t2)
                if m.get_bool('doZeroCon') and mn2 < m.get_int('zeroConPer80'):
                    c2 = 0
                m.set('tMaskTwoGratingContrast', c2)

        # Mask direction only (no con)
        if do_rand_mask_dir and not do_rand_mask_con:
            n = max(1, int(m.get_float('maskOneGratingDirectionStepN', 12)))
            t1 = mn % n
            d1 = m.get_float('tStimOneGratingDirectionDeg') + m.get_float('maskOneGratingDirectionDeg') + m.get_float('maskOneGratingDirectionStepDeg', 30) * t1
            m.set('tMaskOneGratingDirectionDeg', d1)
            if match_s1:
                m.set('tMaskTwoGratingDirectionDeg', d1)
            else:
                n2 = max(1, int(m.get_float('maskTwoGratingDirectionStepN', 12)))
                t2 = mn2 % n2
                d2 = m.get_float('tStimTwoGratingDirectionDeg') + m.get_float('maskTwoGratingDirectionDeg') + m.get_float('maskTwoGratingDirectionStepDeg', 30) * t2
                m.set('tMaskTwoGratingDirectionDeg', d2)

        # Mask phase (no con)
        if do_rand_mask_phase and not do_rand_mask_con:
            n = max(1, int(m.get_float('maskOneGratingPhaseStepN', 12)))
            t1 = mn % n
            m.set('tMaskOneGratingPhaseDeg', m.get_float('maskOneGratingPhaseDeg') + m.get_float('maskOneGratingPhaseStepDeg', 30) * t1)
            if match_s1:
                m.set('tMaskTwoGratingPhaseDeg', m.get_float('tMaskOneGratingPhaseDeg'))
            else:
                n2 = max(1, int(m.get_float('maskTwoGratingPhaseStepN', 12)))
                t2 = mn % n2
                m.set('tMaskTwoGratingPhaseDeg', m.get_float('maskTwoGratingPhaseDeg') + m.get_float('maskTwoGratingPhaseStepDeg', 30) * t2)

        # Mask direction + contrast crossed
        if do_rand_mask_dir and do_rand_mask_con:
            dir_n = max(1, int(m.get_float('maskOneGratingDirectionStepN', 12)))
            con_n = max(1, int(m.get_float('maskOneGratingContrastStepN', 5)))
            t1 = mn % (dir_n * con_n)
            if t1 < dir_n:
                d1 = m.get_float('tStimOneGratingDirectionDeg') + m.get_float('maskOneGratingDirectionDeg') + m.get_float('maskOneGratingDirectionStepDeg', 30) * t1
                c1 = m.get_float('maskOneGratingContrast')
            else:
                d1 = m.get_float('tStimOneGratingDirectionDeg') + m.get_float('maskOneGratingDirectionDeg') + m.get_float('maskOneGratingDirectionStepDeg', 30) * (t1 % dir_n)
                c1 = m.get_float('maskOneGratingContrast') * pow(m.get_float('maskOneGratingContrastStepLog', 2), t1 // dir_n)
            m.set('tMaskOneGratingDirectionDeg', d1)
            m.set('tMaskOneGratingContrast', c1)
            if match_s1:
                m.set('tMaskTwoGratingDirectionDeg', d1)
                m.set('tMaskTwoGratingContrast', c1)

        # Mask phase + contrast crossed
        if do_rand_mask_phase and do_rand_mask_con:
            phase_n = max(1, int(m.get_float('maskOneGratingPhaseStepN', 12)))
            con_n = max(1, int(m.get_float('maskOneGratingContrastStepN', 5)))
            t1 = mn % (phase_n * con_n)
            if t1 < con_n:
                m.set('tMaskOneGratingPhaseDeg', m.get_float('maskOneGratingPhaseDeg'))
            else:
                m.set('tMaskOneGratingPhaseDeg', m.get_float('maskOneGratingPhaseDeg') + m.get_float('maskOneGratingPhaseStepDeg', 30) * (t1 // con_n))
            if match_s1:
                m.set('tMaskTwoGratingPhaseDeg', m.get_float('tMaskOneGratingPhaseDeg'))

    def _apply_matrix_randomization(self, m, sn, sn2, is_b2):
        """Apply matrix (crossed parameter) randomization."""
        match_s1 = m.get_bool('matchStimOneParameters')
        do_rand_dir = m.get_bool('tDoRandDir')
        do_rand_sf = m.get_bool('tDoRandSF')
        do_rand_pos = m.get_bool('tDoRandPos')
        do_rand_isi = m.get_bool('tDoRandISITime')
        do_rand_on = m.get_bool('tDoRandStimOnTime')
        do_rand_con = m.get_bool('tDoRandCon')
        do_mask = m.get_bool('tDoMask')

        # Dir x SF
        if do_rand_dir and do_rand_sf and not do_rand_pos:
            dir_n = max(1, int(m.get_float('stimOneGratingDirectionStepN', 12)))
            sf_n = max(1, int(m.get_float('stimOneGratingSpatialFreqStepN', 5)))
            t1 = sn % (dir_n * sf_n)
            if t1 < dir_n:
                d1 = m.get_float('stimOneGratingDirectionDeg') + m.get_float('stimOneGratingDirectionStepDeg', 30) * t1
                sf1 = m.get_float('stimOneGratingSpatialFreqCPD')
            else:
                d1 = m.get_float('stimOneGratingDirectionDeg') + m.get_float('stimOneGratingDirectionStepDeg', 30) * (t1 % dir_n)
                sf1 = m.get_float('stimOneGratingSpatialFreqCPD') * pow(m.get_float('stimOneGratingSpatialFreqStepLog', 2), t1 // dir_n)
            m.set('tStimOneGratingDirectionDeg', d1)
            m.set('tMaskOneGratingDirectionDeg', d1 + m.get_float('maskOneGratingDirectionDeg'))
            m.set('tStimOneGratingSpatialFreqCPD', sf1)
            if match_s1:
                m.set('tStimTwoGratingDirectionDeg', d1)
                m.set('tStimTwoGratingSpatialFreqCPD', sf1)

        # Dir x ISI
        if do_rand_isi and do_rand_dir:
            dir_n = max(1, int(m.get_float('stimOneGratingDirectionStepN', 12)))
            isi_n = max(1, int(m.get_float('isiTimeN', 5)))
            t1 = sn % (dir_n * isi_n)
            if t1 < dir_n:
                d1 = m.get_float('stimOneGratingDirectionDeg') + m.get_float('stimOneGratingDirectionStepDeg', 30) * t1
                m.set('tISITimeMs', m.get_float('isiTimeMs'))
            else:
                d1 = m.get_float('stimOneGratingDirectionDeg') + m.get_float('stimOneGratingDirectionStepDeg', 30) * (t1 % dir_n)
                m.set('tISITimeMs', m.get_float('isiTimeMs') * pow(m.get_float('isiTimeStepLog', 2), t1 // dir_n))
            m.set('tStimOneGratingDirectionDeg', d1)
            if match_s1:
                m.set('tStimTwoGratingDirectionDeg', d1)

        # Dir x Contrast
        if do_rand_con and do_rand_dir:
            dir_n = max(1, int(m.get_float('stimOneGratingDirectionStepN', 12)))
            con_n = max(1, int(m.get_float('stimOneGratingContrastStepN', 5)))
            t1 = sn % (dir_n * con_n)
            if t1 < dir_n:
                d1 = m.get_float('stimOneGratingDirectionDeg') + m.get_float('stimOneGratingDirectionStepDeg', 30) * t1
                c1 = m.get_float('stimOneGratingContrast')
            else:
                d1 = m.get_float('stimOneGratingDirectionDeg') + m.get_float('stimOneGratingDirectionStepDeg', 30) * (t1 % dir_n)
                c1 = m.get_float('stimOneGratingContrast') * pow(m.get_float('stimOneGratingContrastStepLog', 2), t1 // dir_n)
            m.set('tStimOneGratingDirectionDeg', d1)
            m.set('tStimOneGratingContrast', c1)
            if match_s1:
                m.set('tStimTwoGratingDirectionDeg', d1)
                m.set('tStimTwoGratingContrast', c1)

        # Matrix mask handling
        if do_mask and not m.get_bool('doRandMaskDir') and not m.get_bool('doRandMaskCon'):
            m.set('tMaskOneGratingDirectionDeg', m.get_float('tStimOneGratingDirectionDeg') + m.get_float('maskOneGratingDirectionDeg'))
            m.set('tMaskTwoGratingDirectionDeg', m.get_float('tStimTwoGratingDirectionDeg') + m.get_float('maskTwoGratingDirectionDeg'))
            mult_val = random.random()
            fract = m.get_float('fractMaskTrials')
            if mult_val >= fract:
                m.set('tDoMask', 0)
                m.set('tMaskOneGratingContrast', 0)
                m.set('tMaskTwoGratingContrast', 0)
            else:
                m.set('tDoMask', 1)
                m.set('tMaskOneGratingContrast', m.get_float('tStimOneGratingContrast'))
                m.set('tMaskTwoGratingContrast', 0)
                if m.get_bool('tStimTwoDoVisualStim'):
                    m.set('tMaskTwoGratingContrast', m.get_float('tStimTwoGratingContrast'))

    def _post_randomization_fixups(self, m, sn, sn2, is_b2):
        """Post-randomization fixups matching MWEL lines 2639-2709."""
        # Reset mask phase when stim or mask contrast is 0
        if m.get_float('tMaskOneGratingContrast') == 0 or m.get_float('tStimOneGratingContrast') == 0:
            m.set('tMaskOneGratingPhaseDeg', m.get_float('stimOneGratingPhaseDeg'))
            m.set('tMaskTwoGratingPhaseDeg', m.get_float('stimTwoGratingPhaseDeg'))

        # Default mask TF to stim TF if not randomizing
        if not m.get_bool('doRandMaskTF'):
            m.set('tMaskOneGratingTemporalFreqCPS', m.get_float('stimOneGratingTemporalFreqCPS'))
            m.set('tMaskTwoGratingTemporalFreqCPS', m.get_float('stimTwoGratingTemporalFreqCPS'))

        # Zero contrast when visual stim disabled
        if not m.get_bool('tStimOneDoVisualStim'):
            m.set('tStimOneGratingContrast', 0)
            m.set('tMaskOneGratingContrast', 0)
        if not m.get_bool('tStimTwoDoVisualStim') and not m.get_bool('doTwoStimTogether'):
            m.set('tStimTwoGratingContrast', 0)
            m.set('tMaskTwoGratingContrast', 0)

        # Clamp contrasts to [0, 1]
        for v in ['tStimOneGratingContrast', 'tMaskOneGratingContrast',
                   'tStimTwoGratingContrast', 'tMaskTwoGratingContrast']:
            val = m.get_float(v)
            m.set(v, max(0.0, min(1.0, val)))

        # If stim + mask > 1, zero the mask
        if m.get_float('tStimOneGratingContrast') + m.get_float('tMaskOneGratingContrast') > 1:
            m.set('tMaskOneGratingContrast', 0)
        if m.get_float('tStimTwoGratingContrast') + m.get_float('tMaskTwoGratingContrast') > 1:
            m.set('tMaskTwoGratingContrast', 0)

    # ── Private: end-of-trial selection management ──

    def _do_end_trial(self):
        """End-of-trial actions: accept selections, advance or reset."""
        # svStimNumber
        self.sv_stim.accept()
        self.n_stim_accepted += 1
        if self.n_stim_accepted >= 80:
            self.sv_stim.reset()
            self.n_stim_accepted = 0
        else:
            self.sv_stim.next()

        # svStimNumber2
        self.sv_stim2.accept()
        self.n_stim_accepted2 += 1
        if self.n_stim_accepted2 >= 80:
            self.sv_stim2.reset()
            self.n_stim_accepted2 = 0
        else:
            self.sv_stim2.next()

        # svMaskNumber
        self.sv_mask.accept()
        self.n_mask_accepted += 1
        if self.n_mask_accepted >= 80:
            self.sv_mask.reset()
            self.n_mask_accepted = 0
        else:
            self.sv_mask.next()

        # svMaskNumber2
        self.sv_mask2.accept()
        self.n_mask_accepted2 += 1
        if self.n_mask_accepted2 >= 80:
            self.sv_mask2.reset()
            self.n_mask_accepted2 = 0
        else:
            self.sv_mask2.next()

    def _reset_selections(self):
        """Reset all selection variables."""
        self.sv_stim.reset()
        self.sv_stim2.reset()
        self.sv_mask.reset()
        self.sv_mask2.reset()
        self.n_stim_accepted = 0
        self.n_stim_accepted2 = 0
        self.n_mask_accepted = 0
        self.n_mask_accepted2 = 0
