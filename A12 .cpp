#include <bits/stdc++.h>
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

static const double epsilon1 = 0.00001;
static const double epsilon2 = 0.00000001;

static size_t m;
static size_t n;

struct variable
{
  size_t label;
  double value;
};

struct eta
{
  size_t col;
  vector<double> values;
};
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
void printMatrix(double *matrix, size_t r, size_t c)
{
  for (size_t row = 0; row < r; ++row)
  {
    for (size_t col = 0; col < c; ++col)
    {
      printf("%10.3f ", matrix[row * (m + n) + col]);
    }
    printf("\n");
  }
}

void printLPInfo(double objFuncCoeff[], variable b[], double *matrix)
{
  printf("m = %lu ", m);
  printf("n = %lu \n", n);

  printf("c = ");

  for (size_t i = 0; i < m + n; ++i)
  {
    printf("%10.3f ", objFuncCoeff[i]);
  }

  printf("\nb = ");

  for (size_t i = 0; i < m; ++i)
  {
    printf("%10.3f ", b[i].value);
  }

  printf("\nA = \n");

  printMatrix(matrix, m, n + m);
};

void printVariables(size_t nonbasic[], variable b[])
{
  printf("N = { ");

  for (size_t i = 0; i < n; ++i)
  {
    printf("x%lu ", nonbasic[i] + 1);
  }

  printf("} B = { ");

  for (size_t i = 0; i < m; ++i)
  {
    printf("x%lu ", b[i].label + 1);
  }

  printf("}\n");
};

void printBbar(variable b[])
{
  printf("bbar = ");

  for (size_t i = 0; i < m; ++i)
  {
    printf("%10.3f ", b[i].value);
  }

  printf("\n");
}

void printFinalVariables(variable b[], size_t nonbasic[])
{
  double varValues[m + n];

  for (size_t row = 0; row < m; ++row)
  {
    varValues[b[row].label] = b[row].value;
  }

  for (size_t col = 0; col < n; ++col)
  {
    varValues[nonbasic[col]] = 0.0;
  }

  printf("Decision variables: ");

  for (size_t i = 0; i < n; ++i)
  {
    printf("x%lu = %5.3f ", i + 1, varValues[i]);
  }

  printf("\nSlack variables: ");

  for (size_t i = n; i < m + n; ++i)
  {
    printf("x%lu = %5.3f ", i + 1, varValues[i]);
  }

  printf("\n");
}

void printFamilyOfSolutions(variable b[], size_t nonbasic[], vector<double> d, double largestCoeff, size_t enteringLabel, double z)
{

  struct varInfo
  {
    size_t label;
    double value;
    size_t row_in_basis;
    bool isBasic;
  };

  varInfo info[m + n];

  for (size_t row = 0; row < m; ++row)
  {
    variable v = b[row];
    info[v.label] = {
        v.label,
        v.value,
        row,
        true};
  }

  for (size_t col = 0; col < n; ++col)
  {
    info[nonbasic[col]] = {
        nonbasic[col],
        0.0,
        0,
        false};
  }

  printf("Decision variables: ");

  for (size_t i = 0; i < n; ++i)
  {
    printf("x%lu = %5.3f ", i + 1, info[i].value);
    if (info[i].isBasic)
    {
      printf("+ %5.3fx%lu ", d[info[i].row_in_basis] * (-1.0), enteringLabel + 1);
    }
  }

  printf("\nSlack variables: ");

  for (size_t i = n; i < m + n; ++i)
  {
    printf("x%lu = %5.3f ", i + 1, info[i].value);
    if (info[i].isBasic)
    {
      printf("+ %5.3fx%lu ", d[info[i].row_in_basis] * (-1.0), enteringLabel + 1);
    }
  }

  printf("\nZ = %5.3f + %5.3fx%lu, ", z, largestCoeff, enteringLabel + 1);

  printf("with x%lu >= 0\n", enteringLabel + 1);
}

bool mComparator(variable v1, variable v2)
{
  return v1.value > v2.value;
}

int main(int argc,
         const char *argv[])
{

  cout << "Enter number of equations : ";
  cin >> m;
  cout << "Enter number of variables : ";
  cin >> n;

  double objFuncCoeff[n + m];

  cout << "Enter coefficients of the objective function: ";
  for (size_t col = 0; col < n; ++col)
  {
    cin >> objFuncCoeff[col];
  }

  for (size_t col = n; col < m + n; ++col)
  {
    objFuncCoeff[col] = 0.0;
  }

  double A[m * (n + m)];

  variable b[m];

  size_t nonbasic[n];

  // for (int i = 0; i < m; i++) {
  //         for (int j = 0; j < n; j++) {
  //             cout << "Input matrix coefficient for (" << i + 1 << ", " << j + 1 << ") : ";
  //             cin >> A_mat[i][j];
  //         }
  //     }
  //     cout << endl;
  printf("Enter the coefficients of the constraints\n");

  for (size_t row = 0; row < m; ++row)
  {
    for (size_t col = 0; col <= n; ++col)
    {
      if (col == n)
      {
        cout << "Enter constant term : ";
        double bRow;
        cin >> bRow;
        variable bVar = {
            n + row,
            bRow};
        b[row] = bVar;
      }
      else
      {
        cout << "Input matrix coefficient for (" << row + 1 << ", " << col + 1 << ") : ";
        cin >> A[row * (m + n) + col];
      }
    }
  }

  for (size_t row = 0; row < m; ++row)
  {
    size_t base = (m + n) * row + n;
    for (size_t col = 0; col < m; ++col)
    {
      if (col != row)
      {
        A[base + col] = 0.0;
      }
      else
      {
        A[base + col] = 1.0;
      }
    }
  }

  for (size_t i = 0; i < n; ++i)
  {
    nonbasic[i] = i;
  }

  printLPInfo(objFuncCoeff, b, A);

  printf("\n\n");

  printVariables(nonbasic, b);

  printBbar(b);

  printf("\n");

  for (size_t row = 0; row < m; ++row)
  {
    if (b[row].value < 0.0)
    {
      printf("The given linear program is infeasible, exiting the program.\n");
      return 0;
    }
  }

  size_t counter = 1;

  vector<eta> pivots{};

  double z = 0.0;

  while (true)
  {
    printf("Iteration%lu\n------------\n", counter);

    vector<double> y(m);

    for (size_t row = 0; row < m; ++row)
    {
      variable v = b[row];
      y[row] = objFuncCoeff[v.label];
    }

    for (auto rIter = pivots.crbegin(); rIter != pivots.crend(); ++rIter)
    {
      eta pivot = *rIter;
      size_t colToChange = pivot.col;
      double yOriginal = y[colToChange];

      for (size_t row = 0; row < pivot.values.size(); ++row)
      {
        if (row != colToChange)
        {
          yOriginal -= pivot.values[row] * y[row];
        }
      }

      double yNew = yOriginal / pivot.values[colToChange];
      y[colToChange] = yNew;
    }

    printf("y = ");

    for (auto iter = y.cbegin(); iter != y.cend(); ++iter)
    {
      printf("%10.3f ", *iter);
    }

    printf("\n");

    vector<variable> cnbars;

    size_t enteringLabel = nonbasic[0];
    double largestCoeff = -1.0;

    printf("cnbar: ");

    for (size_t i = 0; i < n; ++i)
    {
      size_t varLabel = nonbasic[i];
      double cni = objFuncCoeff[varLabel];
      double yai = 0.0;

      for (size_t yIndex = 0; yIndex < m; ++yIndex)
      {
        yai += y[yIndex] * A[yIndex * (m + n) + varLabel];
      }

      double cnbar = cni - yai;

      printf("x%lu %5.3f ", varLabel + 1, cnbar);

      if (cnbar > epsilon1)
      {
        variable v = {
            varLabel,
            cnbar};

        cnbars.push_back(v);

        if (cnbar > largestCoeff)
        {
          largestCoeff = cnbar;
          enteringLabel = varLabel;
        }
      }
    }

    sort(cnbars.begin(), cnbars.end(), mComparator);

    printf("\n");

    if (cnbars.size() == 0)
    {
      printf("\nNo entering var. Optimal value of %5.3f has been reached.\n", z);
      cout<<endl;

      printFinalVariables(b, nonbasic);
      return 0;
    }
    else
    {
      printf("Entering variable is x%lu \n", enteringLabel + 1);
    }

    size_t enteringVariable_index = 0;

    vector<double> d(m);

    size_t leavingLabel;
    size_t leavingRow;
    double smallest_t;

    while (true)
    {

      leavingLabel = -1;
      leavingRow = -1;
      smallest_t = -1;

      if (enteringVariable_index > 0)
      {
        printf("\n\nRechoosing entering variable since the diagonal element in the eta column is close to zero.\n");
      }

      if (enteringVariable_index < cnbars.size())
      {
        enteringLabel = cnbars[enteringVariable_index].label;

        if (enteringVariable_index > 0)
        {
          printf("Entering variable is x%lu", enteringLabel + 1);
          cout<<endl;

        }
      }
      else
      {
        printf("\nNo entering var. Optimal value of %5.3f has been reached.", z);
        cout<<endl;
        printFinalVariables(b, nonbasic);
        return 0;
      }
      cout<<endl;

      for (size_t row = 0; row < m; ++row)
      {
        d[row] = A[row * (m + n) + enteringLabel];
      }

      for (auto iter = pivots.cbegin(); iter != pivots.cend(); ++iter)
      {
        eta pivot = *iter;
        size_t rowToChange = pivot.col;
        double dOriginal = d[rowToChange];

        d[rowToChange] = dOriginal / pivot.values[rowToChange];

        for (size_t row = 0; row < d.size(); ++row)
        {
          if (row != rowToChange)
          {
            d[row] = d[row] - pivot.values[row] * d[rowToChange];
          }
        }
      }

      printf("d = ");

      for (auto iter = d.cbegin(); iter != d.cend(); ++iter)
      {
        printf("%5.3f ", *iter);
      }

      printf("\n");

      for (size_t row = 0; row < d.size(); ++row)
      {
        if (d[row] > 0.0)
        {
          leavingLabel = b[row].label;
          leavingRow = row;
          smallest_t = b[row].value / d[row];
        }
      }

      if (leavingLabel == -1)
      {
        printf("\nThe given LP is unbounded. The family of solutions is:\n");
        printFamilyOfSolutions(b, nonbasic, d, largestCoeff, enteringLabel, z);
        return 0;
      }

      printf("ratio: ");

      for (size_t row = 0; row < d.size(); ++row)
      {
        if (d[row] < 0.0)
        {
          continue;
        }

        double t_row = b[row].value / d[row];

        if (t_row >= 0.0)
        {
          printf("x%lu %5.3f ", b[row].label + 1, t_row);
        }

        if (t_row < smallest_t)
        {
          leavingLabel = b[row].label;
          leavingRow = row;
          smallest_t = t_row;
        }
      }

      if (d[leavingRow] > epsilon2)
      {
        printf("\nLeaving variable is x%lu\n", leavingLabel + 1);
        break;
      }
      else
      {
        enteringVariable_index++;
        continue;
      }
    }

    variable enteringVar = {
        enteringLabel,
        smallest_t};
    b[leavingRow] = enteringVar;

    for (size_t row = 0; row < sizeof(b) / sizeof(b[0]); ++row)
    {
      if (row != leavingRow)
      {
        b[row].value -= d[row] * smallest_t;
      }
    }

    eta pivot = {
        leavingRow,
        d};
    pivots.push_back(pivot);

    printf("E%lu = column %lu: ", counter, leavingRow);

    for (auto iter = d.cbegin(); iter != d.cend(); ++iter)
    {
      printf("%5.3f ", *iter);
    }

    printf("\n");

    nonbasic[enteringLabel] = leavingLabel;

    printVariables(nonbasic, b);

    printBbar(b);
    cout<<endl;
    printf("\nCoefficient of entering variable: %5.3f\nAmount increased for the entering variable is: %5.3f\n", largestCoeff, smallest_t);
    double increasedValue = largestCoeff * smallest_t;

    printf("Increased value: %5.3f\n", increasedValue);

    double originalZ = z;

    z += increasedValue;

    printf("Value of the objective function changed from %5.3f to %5.3f\n\n\n", originalZ, z);

    counter++;
  }

  return 0;
}