#include <bits/stdc++.h>
#include<iostream>
using namespace std;

double M = 10000000000.0;
set<int> artificial;
vector<int> sur, art, slack;
vector<int>basVar;
vector<vector<double>>A;
vector<double>objFunc,rhs;
vector<double>z;
int num_var,num_eqn;
int iter = 0;
bool unbounded = false;
bool not_possible = false;
bool degeneracy = false;
bool non_existing = false;
bool optimal = false;
bool alt_soln = false;
bool isMaximization;
bool method_not_applicable = false;

pair<int,int> deg_sol(vector<double>& ratio, double minsol,vector<int>&basv,int key, vector<vector<double>>& A){
    set<int> idxs;
    for(int i=0;i<ratio.size();i++){
        if(ratio[i]==minsol){
            idxs.insert(i);
        }
    }
    for(auto &x:idxs){
        if(artificial.find(x+1)!=artificial.end()){
        double mn = 1/(A[x][key]);
            return {x,mn};
        }
    }
    int fidx = *idxs.rbegin();
    double mn = 1/(A[fidx][key]);
    return {fidx,mn};

}
void simplex(vector<vector<double>>& A, vector<double>& B, vector<double>& C,vector<int>& basis) 
{
    double sum;
    int m = A.size();
    int n = A[0].size();
    // int i, j;

    // vector<double> cb(N, 0.0); 
    // for (int i = 0; i < n; i++){
    //     cb[i] = C[i]; 
    // }
    vector<double> sol(m), ratio(m);  
    int cnt=0;
    int keycol, keyrow;
    double maxincmz = 0.0;  
    double mininsol = INT_MAX; 
    vector<int> basv(m);  
    vector<int> basv_a(m);
    vector<double> sol_a(m);
    vector<double> Z(n, 0),CminusZ(n, 0);  
    vector<double> keyrowval(n, 0), keycolval(m, 0);  
    for (int i = 0; i < m; i++){
        basv[i] = basis[i]-1;
        cout<<basis[i]<<" ";
        // if(artificial.find(basis[i])==artificial.end())sol[i]= 
        sol[i] = B[i];   
    }
    // return;
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;
    
    bool degen = false;
    int entering=-1, leaving=-1;
    while (check)
    {
        if (iter >= 15) 
            break;
        if (iter != 0) 
        {
            basv[keyrow] = keycol; 
            if(iter>0){
                entering = keycol+1;
                cout<<endl<<"Entering variable: "<<entering<<endl;

            }
            for (int i = 0; i < m; i++)
            {
                if (i == keyrow){
                    sol[i] = sol[i]/keyrowval[keycol]; 
                }
                else{
                    sol[i] = sol[i]-(keycolval[i]*solkey)/keyrowval[keycol]; 
                }

                for (int j = 0; j < (n); j++){
                    if (i == keyrow){
                        A[i][j] = A[i][j]/keyrowval[keycol]; 
                    }
                    else {A[i][j] = A[i][j]-(keycolval[i]*keyrowval[j])/keyrowval[keycol];}
                    
                }
            }
        }
        
        check = 0;
        maxincmz = -90000.0;
        for (int i = 0; i < A[0].size(); i++){
            sum = 0;
            for (int j = 0; j < basv.size(); j++){
                sum += C[basv[j]] * A[j][i];
            }
            Z[i] = sum;
            CminusZ[i] = C[i] - Z[i]; 
            if (CminusZ[i] > 0){
                check = 1;
            }
            if (CminusZ[i]>maxincmz){
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (int i = 0; i < m; i++){
            sum += C[basv[i]] * sol[i]; 
        }
        zsol = sum;
        mininsol = INT_MAX;
        for (int i = 0; i < m; i++){
            keycolval[i] = A[i][keycol]; 
        }
        for (int i = 0; i < m; i++){
            if(keycolval[i]<=0){
                ratio[i]=INT_MAX;
                continue;
            }
            ratio[i] = sol[i] / keycolval[i]; 
            if (ratio[i] < mininsol){
                mininsol = ratio[i];
                keyrow = i;
                leaving = basv[i]+1;
                // cout<<i<<" "<<basv[i]+1<<endl;
            }
            else if(ratio[i]==mininsol){
                degen = true;
            }
        }
        // if(degen){
        // cout<<endl<<"Degeneracy Found"<<endl;
        //     auto it = deg_sol(ratio,mininsol,basv,keycol,A);
        //     keyrow=it.first;
        //     mininsol=it.second;
        //     leaving = basv[keyrow]+1;

        //     // cout<<"MIN : "<<keyrow<<endl;
        // }
        for (int i = 0; i < (n); i++){
            keyrowval[i] = A[keyrow][i];
        }
        solkey = sol[keyrow];
        
        cout << "\n\nIteration no: " << iter << "\n";
        cout << "cj's";
        for (int i = 0; i < (n); i++)
        {
            cout << setw(8) << C[i];
        }
        cout<<endl;
        cout<<"BV"<<setw(8)<<"C_B"<<setw(8)<<"X_B"<<setw(8);
        for(int i=0;i<n;i++){
            if(i!=n-1)
                cout<<"X"<<i+1<<setw(8);
            else 
                cout<<"X"<<i+1;
        }
        cout<<setw(8)<<"Ratio";
        cout<<endl;
        for (int i = 0; i < m; i++)
        {
            cout << "x"<<basv[i]+1<<setw(8)<<C[basv[i]]<<setw(8);
            cout << sol[i] << setw(8);
            for (int j = 0; j < (n); j++)
            {
                cout << A[i][j] << setw(8);
            }
            cout<< ratio[i] << "\n";
        }
        cout<<endl;
        cout <<"Z_i"<<setw(8);
        for (int i = 0; i < (n); i++)
        {
            cout << Z[i] << setw(10);
        }
        cout<<endl;
        cout << "\nZ-i - C_i"<<setw(8);
        set<int> idxs;
        set<int> bs;
        for(auto &x:basv)bs.insert(x);
        for (int i = 0; i < (n); i++)
        {
            if(CminusZ[i]==0){
                if(!check&&bs.find(i)==bs.end()){
                    // alt=true;
                    check=1;
                    cnt++;
                    keycol=i;
                    if(cnt==1){
                        basv_a=basv;
                        sol_a=sol;
                    }
                    if(cnt==2)check=0;
                    // return;
                }
            }
            cout << -1*CminusZ[i] << setw(8);
            if(-1*CminusZ[i]<0&&CminusZ[i]!=0)idxs.insert(i);
        }
        for(auto &x:idxs){
            int cnt=0;
            for(int i=0;i<A.size();i++){
                if(A[i][x]<0)cnt++;
            }
            if(cnt==A.size()){
                cout<<endl<<endl<<"Unbounded Solution"<<endl;
                return;
            }
        }
        if(leaving!=entering)cout<<endl<<"Leaving variable: "<<leaving<<endl;
        cout << "\nMinimum ratio is : " << mininsol << " coming at pivot row : " << keyrow + 1 << "\n";
        cout << "Minimum Z-i - C_i is : " << -maxincmz << " coming at pivot column: " << keycol + 1 << "\n";
        // cout<<"BFS : (";
        // for(int i=0;i<n;i++){
        //     if(i!=n-1)cout<<sol[i]<<", ";
        //     else cout<<sol[i]<<")";
        // }
        cout<<endl;
        cout << "Value of z is :" << zsol << "\n";
        iter++;
    }
    // for (int i = 0; i < m; i++)
    // {
    //     if(artificial.find(basv[i]+1)!=artificial.end()){
    //         if(sol[i]>0){
    //             cout<<endl<<"INFEASIBLE SOLUTION"<<endl;
    //             return;
    //         }
    //     }
    // }
    if(degen)
            cout<<endl<<"Degeneracy Found"<<endl;
    // if(alt){
    //     cout<<endl<<"Alternate Solution exists"<<endl;
    //     for (int i = 0; i < m; i++)
    //     {
    //         cout << " x_ " << basv_a[i] + 1 << " = " << sol_a[i] << " ";
    //     }
    // }
    cout<<endl;
    cout << "\n The final optimal values are : ";
    for (int i = 0; i < m; i++)
    {
        cout << " x_ " << basv[i] + 1 << " = " << sol[i] << " ";
    }
    cout << " And rest all are 0\n And the optimal value of Z is : " << zsol << "\n";
    return;
}

void disTable(vector<int>& varRes) {
    int conCount = rhs.size();
    
    cout << endl << "Objective Function: "<< endl << "Maximize: ";
    int ii = 1;
    
    for (int j = 0; j < objFunc.size(); j++) {
        if(objFunc[j]>=0)cout << "+" << objFunc[j] << "x" << ii++ << " ";
        else cout << objFunc[j] << "x" << ii++ << " ";
    }
    cout << endl;
    cout << "Constraints: " << endl;
    for(int i = 0; i < conCount; i++){
        int cnt = 1;
        for(auto j: A[i]){
            if(j==0) j = 0;
            if(j<0)cout << j << "x" << cnt << " ";
            else cout << "+" << j << "x" << cnt << " ";
            cnt++;
        }
        cout << " = "  << rhs[i] << endl;
    }
    
    cout << "And all variables >= 0" << endl << endl;

    cout << "Surplus variables: ";
    for (auto& var : sur) cout << "x" << var << " ";
    cout << endl;

    cout << "Slack variables: ";
    for (auto& var : slack) cout << "x" << var << " ";
    cout << endl;

    cout << "Artificial variables: ";
    for (auto& var : art) cout << "x"<< var << " ";
    cout << endl;

    cout << "Basis variables: ";
    for (auto& var : basVar) cout << "x"<< var << " ";
    cout << endl;
}

void stndForm(vector<int>& consTyp, vector<int>& varRes, int& index) {
    int conCount = rhs.size();
    int iniVar = A[0].size();
    

    for (int i = 0; i < conCount; i++) {
        if(consTyp[i] == 1){
            for(auto &coeff: A[i]) coeff*=-1;
            consTyp[i] = 2;
            rhs[i]*=-1;
        }
       
        slack.push_back(index);
        basVar.push_back(index);
        index++;

        for(int j = 0; j < conCount; j++){
            if(i==j){
                A[j].push_back(1);
            }
            else{
                A[j].push_back(0);
            }
        }
    }

    while(objFunc.size()!=(int)A[0].size()) objFunc.push_back(0);
    num_var = A[0].size();
    disTable(varRes);
}

void printTable(int enter) {
    cout << endl;
    z.assign(num_var, 0.0);

    for (int i = 0; i < num_var; i++) {
        for (int j = 0; j < num_eqn; j++) {
            z[i] += A[j][i] * objFunc[basVar[j] - 1];
        }
        z[i] -= objFunc[i];
    }

    cout << setw(22) << " ";
    for (int i = 0; i < num_var; i++) {
        cout << setw(10) << objFunc[i];
    }
    cout << endl;

    cout << setw(1) << "Xb" << setw(10) << "Cb" << setw(10) << "b";;
    for (int i = 0; i < num_var; i++) {
        cout << setw(9) << "X" << i + 1;
    }
    cout << setw(15) << "min-ratio" << endl;

    for (int i = 0; i < (int)basVar.size(); i++) {
        cout << "x" << basVar[i];
        cout << setw(10) << objFunc[basVar[i] - 1];
        cout << setw(10) << rhs[i];
        for (int j = 0; j < num_var; j++) {
            if(A[i][j]!=0)cout << setw(10) << A[i][j];
            else cout << setw(10) << abs(A[i][j]);
        }
        if(enter!=-1 && A[i][enter]>0) cout << setw(10) << rhs[i]/A[i][enter] << endl;
        else cout << setw(10) << "-" << endl;
    }

    cout << setw(5) << "Xbt.Cb-Cj" << setw(13) << " ";
    for (int i = 0; i < num_var; i++) {
        cout << setw(10) << z[i];
    }
    cout << endl << endl << endl;
}

void print_iter_info(int enter,int leave){
    printTable(enter);
    if(degeneracy){
        degeneracy=false;
        cout << " *** DEGENERACY OCCURED *** " << endl << endl;
    }
    //basic feasible solution
    vector<double>ans(num_var,0.0);
    for(int i = 0; i < num_eqn; i++) ans[basVar[i]-1] = rhs[i];
    cout << "BASIC FEASIBLE SOLUTION: ";
    for(auto i: ans) cout << i << " ";
    cout << endl;

    //basic variables
    cout << "BASIC VARIABLES: ";
    vector<bool>var(num_var,0);
    for(auto it: basVar){
        cout << "x" << it << " ";
        var[it-1]=1;
    }
    cout << endl;

    //non-basic variables
    cout << "NON-BASIC VARIABLES: ";
    for(int i = 0; i < num_var; i++){
        if(var[i]) continue;
        else cout << "x" << i+1 << " ";
    }
    cout << endl;
    


    //incoming variable
    cout << "INCOMING VARIABLE: ";
    if(enter!=-1)cout << "x" << enter+1 << endl;
    else cout << " " << endl;

    //outgoing variable
    cout << "OUTGOING VARIABLE: ";
    if(leave!=-1)cout << "x" << basVar[leave] << endl;
    else cout << " " << endl;

    if(enter!=-1)cout << "PIVOT ELEMENT: (x" << basVar[leave] << ",x" << enter+1 << ")" << endl;
    if(enter!=-1)cout << "PIVOT ELEMENT VALUE:" << A[leave][enter] << endl;
    //minimum ratio

    //objective value
    double val = 0.0;
    for(int i = 0; i < num_var; i++){
        val += (ans[i]*objFunc[i]);
    }
    if(!optimal)cout <<"OBJECTIVE VALUE IS : " << val << endl;
    else{
        cout << endl << "*** OPTIMALITY REACHED ***" << endl;
        cout << endl << "OPTIMAL SOLUTION IS: ";
        for(auto i: ans) cout << i << " ";
        cout << endl;

        cout << endl << "**** OPTIMAL VALUE IS : " << val << " *****" << endl;
        if(isMaximization==false){
            cout << endl << "**** OPTIMAL VALUE OF MINIMISATION PROBLEM IS : " << -1*val << " *****" << endl;
        }
    } 
}

int find_entering_var(int index){
    if(index==-1){
        for(int i = 0; i < num_var; i++){
            double val = 0.0;
            for(int j = 0; j < num_eqn; j++){
                val += A[j][i]*objFunc[basVar[j]-1];
            }
            val-=objFunc[i];
            if(val<0) return -2;
        }
        return -1;
    }
    int ind = -1;
    double maxm = -10000000000.0;
    set<int>zero_pos;
    for(int i = 0; i < num_var; i++){
        double val = 0.0;
        for(int j = 0; j < num_eqn; j++){
            val += A[j][i]*objFunc[basVar[j]-1];
        }
        val-=objFunc[i];
        if(val<0) return -2;
        if(A[index][i]>=0) continue;
        val/=A[index][i];
        if(val>=maxm){
            maxm = val;
            ind = i;
        }
    }
    
    return ind;
}

int find_leaving_var(){
    double minm = 0;
    vector<int>rep_ratio;
    for(int i = 0; i < num_eqn; i++){
        double val = rhs[i];
        if(val<minm){
            minm = val;
            rep_ratio.clear();
        }
        if(val==minm){
            rep_ratio.push_back(i);
        }
    }
    if(rep_ratio.size()==0) return -1;
    int ind = *rep_ratio.begin();
    return ind;
}

void simplex(){
    int ind2 = find_leaving_var();
    
    int ind1 = find_entering_var(ind2);
    if(ind1==-2){
        method_not_applicable = true;
        print_iter_info(-1,-1);
        return;
    }
    if(ind1==-1 && ind2!=-1){
        print_iter_info(-1,-1);
        cout << "*** NON-EXISTING FEASIBLE SOLUTION ***" << endl;
        return;
    }
    if(ind1==-1){
        optimal = true;
        print_iter_info(-1,-1);
        return;
    }
    iter++;
    cout << endl << "ITERATION " << iter << ": " << endl;
    print_iter_info(ind1,ind2);
    for(int i = 0; i < num_var; i++){
        if(i==ind1) continue;
        A[ind2][i] /= A[ind2][ind1];
    }
    rhs[ind2]/=A[ind2][ind1];
    A[ind2][ind1] = 1;

    for(int i = 0; i < num_eqn; i++){
        if(i==ind2) continue;
        for(int j = 0; j < num_var; j++){
            if(j==ind1) continue;
            A[i][j]-=(A[ind2][j]*A[i][ind1]);
        }
        rhs[i]-=(rhs[ind2]*A[i][ind1]);
        A[i][ind1] = 0;
    }
    basVar[ind2] = ind1+1;
    simplex();
    return;
}

signed main() {
    bool simpl = true;
    int numVars;
    cout << "Enter number of equations: ";
    cin >> num_eqn;
    cout << "Enter number of variables: ";
    cin >> numVars;
    int index = 1;
    vector<int> varRes(3 * numVars, 0);

    for (int i = 0; i < numVars; i++) {
        cout << "Is variable " << i + 1 << " unrestricted (1) or restricted (0)? ";
        cin >> varRes[i];
    }
    objFunc.assign(numVars,0.0);
    vector<double>obj;
    cout << "Enter coefficients of the objective function: " << endl;

    cout << "Do you want to maximize or minimize (1 for Maximization, 0 for Minimization) : ";
    cin >> isMaximization;

    for (int i = 0; i < numVars; i++){
        double xx; cin>>xx;
        obj.push_back(xx);
        if(varRes[i]==1) obj.push_back(xx*-1);
    }
    objFunc = obj;
    if (!isMaximization) {
        for (double& coeff : objFunc) coeff *= -1;
    }

    rhs.clear();
    A.clear();
    vector<int> consTyp;


    cout << "Specify constraint type (0 for '=', 1 for '>=', 2 for '<='): " << endl;
    cout<< "Input Format is Coeff(A_i); Type of constraint; Value(B_i)"<<endl;
    int equations = num_eqn;
    for (int i = 0; i < num_eqn; i++) {
        cout << "Enter eqn " << i + 1 << ": ";
        vector<double>temp;
        for (int j = 0; j < numVars; j++) {
            double x; cin>>x;
            temp.push_back(x);
            if(varRes[j]) temp.push_back(-x);
        }
        int constraint;
        double rhs_side;
        cin>>constraint;
        cin>>rhs_side;
        A.push_back(temp);
        rhs.push_back(rhs_side);
        if(constraint)consTyp.push_back(constraint);
        else consTyp.push_back(2);

        if(constraint==0){
            A.push_back(temp);
            rhs.push_back(rhs_side);
            consTyp.push_back(1);
            equations++;
        }
        
    }
    for(int j = 0; j < numVars; j++){
        if(varRes[j]) index++;
        index++;
    }
    cout << endl;
    num_eqn = equations;

    stndForm(consTyp, varRes, index);
    cout << endl;   
    simplex();
    if(method_not_applicable){
        cout << "*** Method Not Applicable ***" << endl;
    }


    cout << "TOTAL ITERATIONS: " << iter << endl;
    return 0;
}