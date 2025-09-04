#include <bits/stdc++.h>
#include <iostream>
using namespace std;
set<int> surplus, artificial, slack;
vector<double> orig_func;
int mx = 0;
bool alt=false;
set<string> st;
int N;
void printRes(vector<vector<double>>& A, vector<double>& b, vector<double>& c, vector<vector<string>>& variables, vector<int>& rstrict,  vector<int>& basis);
pair<int,int> deg_sol(vector<double>& ratio, double minsol,vector<int>&basv,int key, vector<vector<double>>& A);

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
        bool degen = false;
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
        if(degen){
            cout<<"Degeneracy Found"<<endl;
            auto it = deg_sol(ratio,mininsol,basv,keycol,A);
            keyrow=it.first;
            mininsol=it.second;
            leaving = basv[keyrow]+1;

            // cout<<"MIN : "<<keyrow<<endl;
        }
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
                    alt=true;
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
    for (int i = 0; i < m; i++)
    {
        if(artificial.find(basv[i]+1)!=artificial.end()){
            if(sol[i]>0){
                cout<<endl<<"INFEASIBLE SOLUTION"<<endl;
                return;
            }
        }
    }
    if(alt){
        cout<<endl<<"Alternate Solution exists"<<endl;
        for (int i = 0; i < m; i++)
        {
            cout << " x_ " << basv_a[i] + 1 << " = " << sol_a[i] << " ";
        }
    }
    cout<<endl;
    cout << "\n The final optimal values are : ";
    for (int i = 0; i < m; i++)
    {
        cout << " x_ " << basv[i] + 1 << " = " << sol[i] << " ";
    }
    cout << " And rest all are 0\n And the optimal value of Z is : " << zsol << "\n";
    return;
}
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

void StandardConversion(vector<vector<double>>& A, vector<double>& b, vector<double>& c, vector<vector<string>>& variables, vector<int>& basis, vector<int>& cType, vector<int>&rstrict, int ct) {
    int m = A.size();    
    int n = A[0].size(); 

    int idx  = ct;
    for (int i = 0; i < m; i++) {
        if (b[i] < 0) {
            b[i] *= -1;

            if(cType[i]==1)cType[i]=2;
            else if(cType[i]==2)cType[i]=1;
            
            for (int j = 0; j < n; j++) {
                A[i][j]*=-1;
            }
        }
        if (cType[i]==1) {
            string temp = "x";
            temp+=(idx+'0');
            for(int j=0;j<m;j++){
                (j==i?A[j].push_back(-1):A[j].push_back(0));
            variables[j].push_back(temp);
            }
            variables[m].push_back(temp);
            surplus.insert(idx);
            // coeff_check[i].insert(idx);
            idx++;
            temp="x";
            temp+=(idx+'0');
            for(int j=0;j<m;j++){
                (j==i?A[j].push_back(1):A[j].push_back(0));
                variables[j].push_back(temp);
            }
            mx=idx;
            artificial.insert(idx);
            variables[m].push_back(temp);
            basis.push_back(idx);
            // coeff_check[i].insert(idx);
            ++idx;
        } else if (cType[i]==2) {
            string temp = "x";
            temp+=(idx+'0');
            for(int j=0;j<m;j++){
                (j==i?A[j].push_back(1):A[j].push_back(0));
                variables[j].push_back(temp);
            }
            variables[m].push_back(temp);
            slack.insert(idx);
            basis.push_back(idx);
            // coeff_check[i].insert(idx);
            mx=idx;
            ++idx;
        } else{
            string temp = "x";
            temp+=(idx+'0');
            for(int j=0;j<m;j++){
                (j==i?A[j].push_back(1):A[j].push_back(0));
                variables[j].push_back(temp);
            } 
            variables[m].push_back(temp);
            artificial.insert(idx);
            basis.push_back(idx);
            // coeff_check[i].insert(idx);
            mx=idx;
            ++idx;
        }
    }
    if(artificial.size()!=0)cout<<endl<<"We are using Big M method"<<endl;
    else cout<<endl<<"We are using Simplex method"<<endl;
    printRes(A,b,c,variables,rstrict,basis);

}

void printRes(vector<vector<double>>& A, vector<double>& b, vector<double>& c,  vector<vector<string>>& var, vector<int>&rstrict,vector<int>& basis) {
    int m = A.size();
    cout<<"Standard Form: "<<endl;
    cout << "Maximize: ";
    for (int j = 0; j < var[m].size(); j++) {
        if(artificial.find(j+1)!=artificial.end())c[j]=-1e7;
        if(j!=var[m].size()-1){
            cout << c[j] << var[m][j] << (c[j+1]>=0 ? " + " : " ");
        }
        else cout << c[j] << var[m][j] ;
    }
    cout << endl << "Subject to:" << endl;
    for (int i = 0; i < m; i++) {
        int n = A[i].size();
        for(int j=0;j<var[i].size();j++){
            string temp1;
            temp1 = var[i][j];
            if(temp1!="")
            st.insert(temp1); 
        }
        for (int j = 0; j < n; j++) {
            
            if(j!=n-1){
                cout << A[i][j] << var[i][j] << (A[i][j+1]>=0 ? " + " : " ");
            }
            else 
                cout << A[i][j] << var[i][j] ;
            // if(j!=n-1){
            //     cout << A[i][j] << (A[i][j+1]>=0 ? " + " : " ");
            // }
            // else 
            //     cout << A[i][j];
        }
        cout << " = " << b[i] << endl;
    }
    cout<<endl;
    simplex(A,b,c,basis);
    return;
}

int main() {
    cout << "Enter number of eqns: ";
    int m;
    cin >> m;
    cout << "Enter number of variables: ";
    int n;
    cin >> n;
    N=n;
    vector<vector<double>> A_mat(m, vector<double>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << "Input matrix coefficient for (" << i + 1 << ", " << j + 1 << ") : ";
            cin >> A_mat[i][j];
        }
    }
    cout << endl;
    
    vector<int>cType(m,0);
    for(int i=0;i<m;i++){
        cout<<"Enter the type of condition for eqn "<<i+1<<" (0 for =, 1 for >=, 2 for <=) : ";
        cin>>cType[i];
    } 
    cout<<endl;
    
    vector<double> val(m);
    for (int i = 0; i < m; i++) {
        cout << "Enter constant terms for eqn " << i + 1 << ": ";
        cin >> val[i];
    }
    cout << endl;

    cout<<"Enter eqn to maximize or minimize: "<<endl;
    vector<double> func(n+m,0);
    for(int i=0;i<n;i++) {double inp;cin>>inp;func[i]=(inp);}
    bool flag = 1;
    cout<<"You want to maximize or minimize(1 for maximize ; 0 for minimize): ";
    cin>>flag;
    
    if(flag==0)for(auto &x:func)x*=-1;

    vector<int> rstrict(3*n,0);
    for(int i=0;i<n;i++){
        cout<<"Is variable "<<i+1<<" is restricted (enter 0) or unrestricted (enter 1) : ";
        cin>>rstrict[i];
    }
    cout<<endl;
    int ct=n+1;
    vector<vector<string>> variables(m+1);
    map<int,string> mp;
    for(int i=0;i<m+1;i++){
        for(int j=0;j<n;j++){
            string temp = "x";
            temp+=(j+'1');
            variables[i].push_back(temp);

        }
    }
    char ch = 'a';
    for(int i=0;i<n;i++){
        if(rstrict[i]){
            string tt = "x";
            tt+=(ct+'0');
            variables[m].push_back(tt);
            func.push_back(0);
            for(int j=0;j<m;j++){
                A_mat[j].push_back(-1*A_mat[j][i]);
                variables[j].push_back(tt);
            }
            ct++;
        }
    }
    orig_func=func;
    vector<int> basis;
    StandardConversion(A_mat,val,func,variables,basis,cType,rstrict,ct);

    return 0;
}