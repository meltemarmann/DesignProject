/***************
 * OPL 22.1.0.0 Model
 * Author: MeltemArman
 * Creation Date: 18 Þub 2023 at 15:54:21
 ***************/
range T = 0..8760;

//Parametreler
float D[T] = ...; //t. saatteki talep
float G_b = 1; // Bataryanin maksimum dolum ve boþaltim hizi (c rate).
float K_pv = 800000; // Batarya kapasitesinin yatirim maliyeti ($/MWh).
float K_w = 1300000; //Gunes enerjisi santralinin yatirim maliyeti ($/MW).
float K_b = 500000; //Ruzgar enerjisi santralinin yatirim maliyeti ($/MW).
float S_min = 0.10; //Bataryada minimum depolanabilecek guc miktari (%).
float S_max = 0.90; //Bataryada maximum depolanabilecek guc miktari (%).
float g_pv[T] = ...; //Gunes enerjisi santralinin t. saatteki kapasite faktoru (%).
float g_w[T] = ...; //Ruzgar enerjisi santralinin t. saatteki kapasite faktoru (%).
float n = 0.80; //Bataryanin dolum ve bosaltim verimlilik orani (%).

//Karar Degiskenleri
dvar float+ c_b; //Bataryanin toplam kapasitesi.
dvar float+ M_pv; //Gunes enerjisi uretim miktarinin carpim katsayisi.
dvar float+ M_w; //Ruzgar enerjisi uretim miktarinin carpim katsayisi.
dvar float+ X_ch[T]; //t. saatte bataryaya doldurulan guc miktari.
dvar float+ X_dis[T]; //t. saatte bataryadan bosaltilan guc miktari.
dvar float+ S[T]; //t. saatte bataryada bulunan guc miktari.
dvar boolean b_ch[T]; //t. saatte batarya dolduruluyorsa 1, diger durumlarda 0 degerini alir.
dvar boolean b_dis[T]; //t. saatte batarya bosaltiliyorsa 1, diger durumlarda 0 degerini alir.
//dvar float+ T_pv; //Herhangi bir saatte gunes enerjisi santralinden uretilen maximum guc miktari.
//dvar float+ T_w; // Herhangi bir saatte ruzgar enerjisi santralinden uretilen maximum guc miktari.

//Amac Fonksiyonu
minimize K_pv*M_pv + K_w*M_w + K_b*c_b;

//Kisitlar
subject to{

forall(t in T){
	g_pv[t]*M_pv + g_w[t]*M_w + X_dis[t]*n  >= D[t] + X_ch[t]/n; //Akis dengeleme ve talep karsilama kisiti.
}
forall(t in T){
	S_min*c_b <= S[t]; //Batarya minimum bosaltilabilecegi miktardan daha fazla bosaltýlamaz.
	S_max*c_b >= S[t]; //Batarya maximum doldurulabilecegi miktardan daha fazla doldurulamaz.
}
forall(t in T){
	X_ch[t] <= G_b*c_b; //Batarya doldurma/bosaltma hizindan daha hizli doldurulamaz.
	X_dis[t] <= G_b*c_b; //Batarya doldurma/bosaltma hizindan daha hizli bosaltilamaz.
}
forall(t in T: t!=0){
	S[t] == S[t-1] + X_ch[t] - X_dis[t]; //Bataryanin doluluk miktarini guncelleme kisiti.
}
forall(t in T){
	b_ch[t] + b_dis[t] <= 1; //Batarya ayni anda doldurulup bosaltilamaz.
}
forall(t in T){ //Binary degiskenleri iliskilendirme kisitlari.
	X_ch[t] <= b_ch[t]*(sum(i in T)D[i]); 
	X_dis[t] <= b_dis[t]*(sum(i in T)D[i]);
}



}
