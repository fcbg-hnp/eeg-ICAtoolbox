# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\icatoolbox_v01.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtWidgets
from toolbox import Ui_MainWindow
import mne
import mne.channels
from read import read_sef
from util import xyz_to_montage, plot_overlay


class MainWindow(Ui_MainWindow):

    def init_variables(self):
        """Initialize variable at start"""
        self.set_tooltips()
        self.pushButton_openxz.setEnabled(False)
        self.lineEdit_xyz.setEnabled(False)
        self.raw = None
        self.n_channels = None
        self.montage = mne.channels.read_montage("standard_1005")
        self.groupBox_advancedparameters.setChecked(False)
        self.collapse()
        self.Openfile_eeg = QtWidgets.QFileDialog(caption='Open file')
        self.Openfile_xyz = QtWidgets.QFileDialog(caption='Open file')
        self.messagebox = QtWidgets.QMessageBox()
        self.reset_variables()
        return()

    def set_tooltips(self):
        self.pushButton_compute.setToolTip('Compute ICA.')
        self.pushButton_apply.setToolTip("Apply current ICA to EEG data")
        self.pushButton_save.setToolTip("Save processed data")
        self.pushButton_plotsources.setToolTip("Plot independent components time serie.\n Click on a source to exclude component from signal reconstruction.\n Click on source label to show component topomap")
        self.pushButton_plottopomaps.setToolTip("""Plot independent components topomaps. \n Click on label to exclude component from reconstruction.\n Click on topomap to show additional informations.""")
        self.pushButton_plotoverlay.setToolTip("Plot Global field power of original EEG signals and reconstructed signals considered current excluded components.")
        self.label_method.setToolTip("Method used to compute ICA solution")
        self.label_maxiter.setToolTip("Maximum number of iterations during ICA fitting. \n If method does not converge after max_iter iteration, it stops and returns current solution. \n If method converge before max_iter, it stops and return the completed solution")
        self.label_randomseed.setToolTip("Random state to initialize ICA estimation for reproducible results.")
        self.label_ncomponents.setToolTip("Controls the number of PCA components from the pre-ICA PCA entering the ICA decomposition.\n Must be < to max_pca_components")
        self.label_maxpcacomponents.setToolTip("The number of components returned by the pre-ICA PCA decomposition.\n Must be < to the number of channels")
        self.label_npcacomponents.setToolTip("The number of PCA components used for re-projecting the decomposed data into sensor space. \n Has to be >= n_components and <= max_pca_components.\n If greater than n_components_, the next n_pca_components minus n_components PCA components \n will be added before restoring the sensor space data.")
        return()

    def reset_variables(self):
        """Reset variable to default """
        self.clean_raw = None
        self.ica = None
        self.save_name = None
        self.fsave = None
        self.validfile = False
        self.validparameters = False
        self.computed = False
        self.applied = False
        self.saved = False
        self.update_compute_label()
        self.update_apply_label()
        self.update_save_label()
        self.activate_groupbox_compute()
        self.activate_compute_button()
        self.activate_groupbox_plot()
        self.activate_groupbox_apply()
        self.activate_groupbox_save()
        self.valid_inputs()
        return()

    def connect_events(self):
        """Coonect functions to the GUI"""
        self.pushButton_openfile.clicked.connect(self.open_eeg_file)
        self.pushButton_openxz.clicked.connect(self.open_xyz_file)
        self.comboBox_montage.activated[str].connect(self.set_montage_from_combobox)
        self.groupBox_advancedparameters.toggled.connect(self.collapse)
        self.groupBox_advancedparameters.toggled.connect(self.set_default_advanced_parameters)
        self.comboBox_methods.activated[str].connect(self.set_method)
        self.spinBox_ncomponents.editingFinished.connect(self.set_ncomponents)
        self.spinBox_seed.editingFinished.connect(self.set_seed)
        self.spinBox_maxiter.editingFinished.connect(self.set_maxiter)
        self.spinBox_maxpcacomponents.editingFinished.connect(self.set_maxpcacomponents)
        self.spinBox_npcacomponents.editingFinished.connect(self.set_npcacomponents)
        self.pushButton_compute.clicked.connect(self.compute)
        self.pushButton_plotsources.clicked.connect(self.plot_sources)
        self.pushButton_plottopomaps.clicked.connect(self.plot_topomaps)
        self.pushButton_plotoverlay.clicked.connect(self.plot_overlay)
        self.pushButton_apply.clicked.connect(self.apply)
        self.pushButton_save.clicked.connect(self.save)
        return()

    def init_parameters(self):
        """Initialise ICA parameters baed on Raw file infos"""
        self.n_channels = len(mne.pick_types(self.raw.info, meg=True, eeg=True))
        self.method = "fastica"
        self.ncomponents = self.n_channels
        self.spinBox_ncomponents.setValue(self.ncomponents)
        self.seed = 42
        self.spinBox_seed.setValue(42)
        self.maxiter = 500
        self.spinBox_maxiter.setValue(500)
        self.maxpcacomponents = self.n_channels
        self.spinBox_maxpcacomponents.setValue(self.maxpcacomponents)
        self.npcacomponents = self.n_channels
        self.spinBox_npcacomponents.setValue(self.npcacomponents)
        self.check_parameters()
        return()

    def collapse(self):
        """Expand/collapse advanced ica parameters menu"""
        if self.groupBox_advancedparameters.isChecked() is True:
            self.groupBox_advancedparameters.setFixedSize(self.groupBox_advancedparameters.sizeHint())
            self.groupBox_advancedparameters.resize(self.groupBox_advancedparameters.sizeHint())
        else:
            self.groupBox_advancedparameters.setFixedSize(self.groupBox_advancedparameters.sizeHint().width(), 20)
            self.groupBox_advancedparameters.resize(self.groupBox_advancedparameters.sizeHint().width(), 20)
        return()

    # Open file
    def open_eeg(self):
        """Open EEG file (support -raw.fif and .sef)"""
        try:
            self.n_channels = None
            self.Raw = None
            filter = "Raw fif(*-raw.fif);;Raw sef (*.sef)"
            self.fname_eeg, self.ext_eeg = self.Openfile_eeg.getOpenFileName(caption='Open file', filter=filter)
            if self.ext_eeg == "Raw fif(*-raw.fif)":
                self.raw = mne.io.read_raw_fif(self.fname_eeg, preload=True)
            elif self.ext_eeg == "Raw sef (*.sef)":
                self.raw = read_sef(self.fname_eeg)
            self.lineEdit_eegfile.setText(self.fname_eeg)
        except Exception as e:
            self.messagebox.setText(str(e))
            self.messagebox.exec()
        return()

    def open_xyz(self):
        """Open electrode postition file (support .xyz)"""
        try:
            self.Openfile_xyz = QtWidgets.QFileDialog(caption='Open file')
            self.fname_xyz, self.ext_xyz = self.Openfile_xyz.getOpenFileName(caption='Open file', filter="*.xyz")
            print(self.ext_xyz)
            if self.ext_xyz == "*.xyz":
                self.set_montage_from_file()
            self.lineEdit_xyz.setText(self.fname_xyz)
            self.valid_inputs()
        except Exception as e:
            self.messagebox.setText("Unable to read montage from file because of error: " + str(e))
            self.messagebox.exec()
        return()

    def set_montage_from_file(self):
        """Set montage to the one contained in the file"""
        self.reset_variables()
        self.montage = xyz_to_montage(self.fname_xyz)
        QtWidgets.QApplication.processEvents()
        self.valid_inputs()
        return()

    def set_montage_from_combobox(self, montage):
        """ Set montage from combobox list"""
        self.reset_variables()
        if str(montage) != "From file":
            self.pushButton_openxz.setEnabled(False)
            self.lineEdit_xyz.setEnabled(False)
            self.montage = mne.channels.read_montage(montage)
            self.lineEdit_xyz.setText("")
            self.fname_xyz = None
            self.ext_xyz = None
        else:
            self.pushButton_openxz.setEnabled(True)
            self.lineEdit_xyz.setEnabled(True)
            self.montage = None
            self.lineEdit_xyz.setText("")
            self.fname_xyz = None
            self.ext_xyz = None
        self.valid_inputs()
        return()

    def valid_inputs(self):
        """Check if both raw and motage are correctly set"""
        if self.raw is not None and self.montage is not None:
            self.validfile = True
            self.init_parameters()
        else:
            self.validfile = False
        self.activate_groupbox_compute()
        return()

    def open_eeg_file(self):
        """Open eeg file and check if inputs are correct"""
        self.reset_variables()
        self.open_eeg()
        self.valid_inputs()
        return()

    def open_xyz_file(self):
        """Open montage file and check if inputs are correct"""
        self.reset_variables()
        self.open_xyz()
        self.valid_inputs()
        self.activate_groupbox_compute()
        return()

    def activate_groupbox_compute(self):
        """Enable compute groupbox if inputs are corrects"""
        if self.validfile is True:
            self.groupBox_setparameters.setEnabled(True)
            self.init_parameters()
        else:
            self.groupBox_setparameters.setEnabled(False)
        return()

    # Parameters
    def check_parameters(self):
        """Checks if ica parameters are correctly set"""
        if self.ncomponents > self.n_channels:
            self.validparameters = False
        elif self.npcacomponents < self.ncomponents:
            self.validparameters = False
        elif self.npcacomponents > self.maxpcacomponents:
            self.validparameters = False
        elif self.maxpcacomponents > self.n_channels:
            self.validparameters = False
        else:
            self.validparameters = True
        self.activate_compute_button()
        return()

    def activate_compute_button(self):
        """Enable compute button if validparameters = True"""
        if self.validparameters is True:
            self.pushButton_compute.setEnabled(True)
        else:
            self.pushButton_compute.setEnabled(False)
        return()

    def set_default_advanced_parameters(self):
        """Set default ica advanced parameters based on Raw info"""
        self.seed = 42
        self.spinBox_seed.setValue(42)
        self.maxiter = 500
        self.spinBox_maxiter.setValue(500)
        self.maxpcacomponents = self.n_channels
        self.spinBox_maxpcacomponents.setValue(self.maxpcacomponents)
        self.npcacomponents = self.n_channels
        self.spinBox_npcacomponents.setValue(self.npcacomponents)
        self.check_parameters()
        return()

    def set_method(self, method):
        """Select ica methods"""
        self.method = method
        return()

    def set_ncomponents(self):
        """Set ncomponents"""
        self.ncomponents = int(self.spinBox_ncomponents.text())
        self.check_parameters()
        return()

    def set_seed(self):
        """Set ica random seed"""
        self.seed = int(self.spinBox_seed.text())
        self.check_parameters()
        return()

    def set_maxiter(self):
        """Set maximum ica iterations"""
        self.maxiter = int(self.spinBox_maxiter.text())
        self.check_parameters()
        return()

    def set_maxpcacomponents(self):
        """Set maximum pca components"""
        self.maxpcacomponents = int(self.spinBox_maxpcacomponents.text())
        self.check_parameters()
        return()

    def set_npcacomponents(self):
        """Set the number of pca components used for ica reconstruction"""
        self.npcacomponents = int(self.spinBox_npcacomponents.text())
        self.check_parameters()
        return()

    def activate_groupbox_plot(self):
        """Enable groubox plot if ica computation is done"""
        if self.computed is True:
            self.groupBox_plot.setEnabled(True)
        else:
            self.groupBox_plot.setEnabled(False)
        return()

    def activate_groupbox_apply(self):
        """Enable groubox apply if ica computation is done"""
        if self.computed is True:
            self.groupBox_apply.setEnabled(True)
        else:
            self.groupBox_apply.setEnabled(False)
        return()

    def compute_ica(self):
        """Compute ica"""
        try:
            raw = self.raw
            raw.set_montage(self.montage)
            ica = mne.preprocessing.ICA(n_components=self.ncomponents,
                                        method=self.method,
                                        random_state=self.seed,
                                        max_iter=self.maxiter,
                                        max_pca_components=self.maxpcacomponents,
                                        n_pca_components=self.npcacomponents)
            picks = mne.pick_types(raw.info, meg=True, eeg=True,
                                   eog=True, exclude='bads')
            print(self.ncomponents, self.method, self.seed, self.maxiter,
                  self.maxpcacomponents, self.npcacomponents, self.montage, self.raw)
            ica.fit(raw, picks=picks, decim=1)
            self.ica = ica
            self.computed = True
        except Exception as e:
            self.messagebox.setText("Unable to compute ica because of error: " + str(e))
            self.messagebox.exec()
            self.computed = False
        return()

    def update_compute_label(self):
        """Update computation info label"""
        if self.computed is True:
            self.label_compute.setText("Done")
        else:
            self.label_compute.setText("None")
        return()

    def compute(self):
        """Compute ica and update GUI"""
        self.applied = False
        self.saved = False
        self.clean_raw = None
        self.update_apply_label()
        self.update_save_label()
        self.activate_groupbox_plot()
        self.activate_groupbox_apply()
        self.activate_groupbox_save()
        self.label_compute.setText("Computing... This may take a while")
        QtWidgets.QApplication.processEvents()
        self.compute_ica()
        self.update_compute_label()
        self.activate_groupbox_plot()
        self.activate_groupbox_apply()
        return()

    def plot_sources(self):
        """Plot sources"""
        try:
            self.ica.plot_sources(self.raw, block=True)
        except Exception as e:
            self.messagebox.setText("Unable to plot sources because of error: " + str(e))
            self.messagebox.exec()
        return()

    def plot_topomaps(self):
        """Plot topomaps"""
        try:
            self.ica.plot_components(inst=self.raw)
        except Exception as e:
            self.messagebox.setText("Unable to plot components because of error: " + str(e))
            self.messagebox.exec()
        return()

    def plot_overlay(self):
        "Plot overlay"
        try:
            plot_overlay(self.raw, self.ica)
        except Exception as e:
            self.messagebox.setText("Unable to plot overlay because of error: " + str(e))
            self.messagebox.exec()
        return()

    def apply_ica(self):
        """Apply ica"""
        raw_copy = self.raw.copy()
        self.ica.apply(raw_copy)
        self.clean_raw = raw_copy
        return()

    def update_apply_label(self):
        """Update apply label info"""
        if self.applied is True:
            self.label_apply.setText("Done   " + str(self.nout) + " components rejected")
        else:
            self.label_apply.setText("None")
        return()

    def apply(self):
        """Apply ica and update GUI"""
        self.label_apply.setText("Computing... This may take a while")
        self.nout = len(self.ica.exclude)
        QtWidgets.QApplication.processEvents()
        self.apply_ica()
        self.applied = True
        self.update_apply_label()
        self.activate_groupbox_save()

    def activate_groupbox_save(self):
        """Enable save groupbox if ica is applied"""
        if self.applied is True:
            self.groupBox_save.setEnabled(True)
        else:
            self.groupBox_save.setEnabled(False)
        return()

    def update_save_label(self):
        """Update save label info"""
        if self.saved is True:
            self.label_save.setText("Saved. " + str(self.fsave))
        else:
            self.label_save.setText("None")
        QtWidgets.QApplication.processEvents()
        return()

    def generate_fname(self):
        """Generate default saving file name based o input file"""
        self.save_name = self.fname_eeg[0:-len(self.ext_eeg)]

    def save(self):
        """Save file"""
        self.generate_fname()
        try:
            self.save_name = self.Openfile_eeg.getSaveFileName(directory=(self.save_name + "_ica-raw.fif"), filter="*-raw.fif")[0]
            self.clean_raw.save(self.save_name, overwrite=True)
            self.saved = True
            self.label_save.setText("Saved. " + str(self.save_name))
            QtWidgets.QApplication.processEvents()
        except Exception as e:
            self.label_save.setText(str(e))
            self.saved = False
            self.label_save.setText(str(e))
            QtWidgets.QApplication.processEvents()

        return()
