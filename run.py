import sys
from PyQt5.QtWidgets import QApplication, QMainWindow
from icatoolbox import MainWindow

if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = MainWindow()
    w = QMainWindow()
    ex.setupUi(w)
    ex.init_variables()
    ex.connect_events()
    w.show()
    sys.exit(app.exec_())
