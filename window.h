#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QHBoxLayout>
#include <QTableWidget>
#include <QtCharts>
#include <QLineEdit>

#include <QDebug>

#include <eigen3/Eigen/Dense>

using namespace Eigen;

class window : public QWidget{
		Q_OBJECT
		public:
				explicit window(QWidget *parent=0);
		private:
				QHBoxLayout *mainLayout;
				QTabWidget *tabWidget;
};

class phaseSpace : public QWidget{
		Q_OBJECT
		public:
				explicit phaseSpace(QWidget *parent=0);
				private slots:
						void updateChart();
				public slots:
						void processVector(VectorXf y);
		private:
				void createChart();
				void createInput();

				QChartView *chartView;
				QChart *phaseChart;
				QValueAxis *axisX, *axisY;
				QScatterSeries *phaseSeries;
				QLineEdit *Ninput, *vbinput, *Linput, *Jinput, *dtinput, *tmaxinput;
				QGroupBox *inputBox;
				QPushButton *goButton;
};

class worker : public QObject{
		Q_OBJECT
		public:
				worker(float tmax, float dr, int N, int J, float L, float vb);
		public slots:
				void passEval();
		signals:
				void returnVector(VectorXf y);	
				void finished();
		private:
				int N, K, J;
				float tmax, dt, L, vb;	
};
#endif // WINDOW_H
