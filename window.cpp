#include "window.h"
#include "plasma.h"

window::window(QWidget *parent):
    QWidget(parent){
	qRegisterMetaType<VectorXf>("VectorXf");
    setWindowTitle("Plasma");
    mainLayout = new QHBoxLayout();
    mainLayout->addWidget(new phaseSpace());

    setLayout(mainLayout);

    resize(900, 600);
    setMinimumHeight(600);
    setMinimumWidth(900);
}

phaseSpace::phaseSpace(QWidget *parent) : QWidget(parent){
    createChart();
    createInput();

    QHBoxLayout *mainLayout = new QHBoxLayout();
    mainLayout->addWidget(chartView);
    mainLayout->addWidget(inputBox);

    setLayout(mainLayout);
}

void phaseSpace::createChart(){
    phaseChart = new QChart();
    phaseChart->legend()->setVisible(false);

    chartView = new QChartView(phaseChart);
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->setMinimumWidth(500);
}

void phaseSpace::createInput(){
    Ninput = new QLineEdit;
    Ninput->setPlaceholderText("N");
    vbinput = new QLineEdit;
    vbinput->setPlaceholderText("vb");
    Linput = new QLineEdit;
    Linput->setPlaceholderText("L");
    Jinput = new QLineEdit;
    Jinput->setPlaceholderText("J");
    dtinput = new QLineEdit;
    dtinput->setPlaceholderText("dt");
    tmaxinput = new QLineEdit;
    tmaxinput->setPlaceholderText("tmax");

    QLabel *nlabel = new QLabel("Number of Electrons");
    QLabel *vblabel = new QLabel("Beam Velocity");
    QLabel *Llabel = new QLabel("Domain of Solution /Debye Lengths");
    QLabel *Jlabel = new QLabel("Number of Grid Points");
    QLabel *dtlabel = new QLabel("Time step \n/Inverse Plasma Frequency");
    QLabel *tmaxlabel = new QLabel("Time");

    QLabel *description = new QLabel("This program simulates the behavior of a 1D plasma beam over a period of time etc.");
    description->setWordWrap(true);

    goButton = new QPushButton("Go!");

    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->addWidget(nlabel,0,0);
    mainLayout->addWidget(Ninput,0,1);
    mainLayout->addWidget(vblabel,1,0);
    mainLayout->addWidget(vbinput,1,1);
    mainLayout->addWidget(Llabel,2,0);
    mainLayout->addWidget(Linput,2,1);
    mainLayout->addWidget(Jlabel,3,0);
    mainLayout->addWidget(Jinput,3,1);
    mainLayout->addWidget(dtlabel,4,0);
    mainLayout->addWidget(dtinput,4,1);
    mainLayout->addWidget(tmaxlabel,5,0);
    mainLayout->addWidget(tmaxinput,5,1);
    mainLayout->addWidget(description,6,0,1,2);
    mainLayout->addWidget(goButton,8,0,1,2,Qt::AlignBottom);

    connect(goButton, SIGNAL(clicked(bool)), this, SLOT(updateChart()));

    inputBox = new QGroupBox("Input");
    inputBox->setLayout(mainLayout);
    inputBox->setFixedWidth(300);
}

void phaseSpace::updateChart(){
		goButton->setEnabled(false);
    int N = Ninput->text().toInt();
    float L = Linput->text().toFloat();
    int J = Jinput->text().toInt();
    float vb = vbinput->text().toFloat();
    float dt = dtinput->text().toFloat();
    float tmax = tmaxinput->text().toFloat();

    if(N == 0 || L == 0 || J==0 || dt == 0)
        return;

    static int flag = 0;
    if(flag == 1){
        delete phaseSeries;
    }
    flag = 1;
	
	auto thread = new QThread;
	auto worker = new class worker(tmax, dt, N, J, L, vb);
	worker->moveToThread(thread);
	connect(thread, SIGNAL(started()), worker, SLOT(passEval()));
	connect(worker, SIGNAL(returnVector(VectorXf)), this, SLOT(processVector(VectorXf)));
	connect(worker, SIGNAL(finished()), thread, SLOT(quit()));
	connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
	connect(worker, SIGNAL(finished()), thread, SLOT(deleteLater()));
	thread->start();	
}

void phaseSpace::processVector(VectorXf y){
		int N = y.size()/2;
	VectorXf v(N), r(N);
    v = y.tail(N);
    r = y.head(N);

    phaseSeries = new QScatterSeries;
    phaseSeries->setUseOpenGL(true);
    phaseSeries->setColor(QColor("black"));
    phaseSeries->setMarkerSize(2);

    for(int i = 0; i < N; i++){
        phaseSeries->append(r(i), v(i));
        //qDebug() <<  "r(i) = " << r(i) << "v(i) =" << v(i);
    }
    phaseChart->addSeries(phaseSeries);
//    phaseSeries->attachAxis(axisX);
//    phaseSeries->attachAxis(axisY);
    phaseChart->createDefaultAxes();
	goButton->setEnabled(true);
}
worker::worker(float tmax, float dt, int N, int J, float L, float vb)
		: tmax(tmax), dt(dt), N(N), J(J), L(L), vb(vb) {}

void worker::passEval(){
    plasma plasma(L,J,N);
    plasma.rDist();
    plasma.vDist(vb);
	VectorXf y(2*N);
    y = plasma.eval(tmax, dt);
	emit returnVector(y);
	emit finished();
}
