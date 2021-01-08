/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <cerrno>	/* Need for Using "ERANGE" */
#include <cmath>	/* Need for Using "HUGE_VAL" */
#include <iostream>
#include <climits>

#include "delfem2/evalmathexp.h"

namespace dfm2 = delfem2;

class COperand : public delfem2::evalmathexp::CCmd
{
public:
  ~COperand() override = default;
  explicit COperand(const double val){ m_Value = val; }
  
  bool DoOperation(std::vector<double>& stack) override{
    stack.push_back(m_Value);
    return true;
  }
  static double GetValue(const int iopr){
    switch(iopr){
      case 0:
        return 3.14159265358979;
      default:
        assert(0);
        return 0.0;
    }
    return 0.0;
  }
  static int GetMaxOprInd(){ return 0; }
  void SetValue(const double& val) override { m_Value = val; }
  static int GetOprInd(const std::string& str1){
    if( str1 == "PI" ){ return 0; }
    else{ return -1; }
    return -1;
  }
private:
  double m_Value;
};

class CBinaryOperator : public delfem2::evalmathexp::CCmd
{
public:
  ~CBinaryOperator() override = default;
  explicit CBinaryOperator(const unsigned int iopr){
    assert( iopr<5 );
    m_iOpr = iopr;
  }
  
  bool DoOperation(std::vector<double>& stack) override{
    const double dright = stack.back();
    stack.pop_back();
    double& dleft = stack.back();
    switch(m_iOpr){
      case 0:
        dleft = dleft+dright; return true;
      case 1:
        dleft = dright-dleft; return true;
      case 2:
        dleft = dleft*dright; return true;
      case 3:
        if( fabs(dleft) < 1.0e-20 ){ dleft=0.0; return false; }
        dleft = dright/dleft; return true;
      case 4:
        dleft = pow(dright,dleft); return true;
      default:
        std::cout << "Error!-->Illegal iOpr " << m_iOpr << std::endl;
        assert(0);
    }
    return false;
  }
  static int MaxOprInd(){ return 4; }
  void SetValue(const double& val) override{
    assert(0);
  };
  unsigned int m_iOpr;
};

class CUnaryOperator : public delfem2::evalmathexp::CCmd
{
public:
  ~CUnaryOperator() override = default;
  explicit CUnaryOperator(const int iopr){
    assert( iopr>=0 && iopr<=MaxOprInd() );
    m_iOpr = iopr;
  }
  bool DoOperation(std::vector<double>& stack) override{
    const double dright = stack.back();
    double& dleft = stack.back();
    switch(m_iOpr){
      case 0: dleft = dright; return true;
      case 1: dleft = -dright; return true;
      case 2: dleft = fabs(dright); return true;
      case 8: dleft = floor(dright); return true;
      case 3: dleft = exp(dright);  return true;
      case 4: dleft = sin(dright);  return true;
      case 5: dleft = cos(dright);  return true;
      case 6:
        if( dright > 1.0e-30 ){ dleft = sqrt(dright); return true; }
        else{ dleft = 0.0; return false; }
      case 7:
        if( dright > 1.0e-30 ){ dleft = log(dright); return true; }
        else{ dleft = 0.0; return false; }
      default:
        std::cout << "Error!-->Illegal iOpr " << m_iOpr << std::endl;
        assert(0);
    }
    return false;
  }
  static int MaxOprInd(){ return 8; }
  static int GetOprInd(const std::string& str1){
    if( str1 == "+" )      return 0;
    else if( str1 == "-" )    return 1;
    else if( str1 == "abs" )  return 2;
    else if( str1 == "floor" )  return 8;
    else if( str1 == "exp" )  return 3;
    else if( str1 == "sin" )  return 4;
    else if( str1 == "cos" )  return 5;
    else if( str1 == "sqrt" )  return 6;
    else if( str1 == "log" )  return 7;
    return -1;
  }
  void SetValue(const double& val) override{
    assert(0);
  }
private:
  int m_iOpr;
};

// -------------------------------------

static void RemoveExpressionSpaceNewline( std::string& exp )
{
  const std::string stemp(exp);
  unsigned int iexp=0;
  for(char itemp : stemp){
    if( itemp == ' ' || itemp == 0x0d || itemp == 0x0a )  continue;
    exp[iexp] = itemp;
    iexp++;
  }
  exp.resize(iexp);
}

static void RemoveExpressionBracket(std::string& exp)
{
  if( exp[0] != '(' || exp[exp.size()-1] != ')' ) return;
  {
    int iBracketDepth = 1;
    for ( size_t ipos=1; ipos<exp.size()-1; ipos++ )
    {
      if ( exp[ipos] == '(' ){ iBracketDepth++; }
      if ( exp[ipos] == ')' ){ iBracketDepth--; }
      if ( iBracketDepth == 0 ) return;
    }
  }
  {
    std::string stemp(exp);
    for ( size_t ipos=1; ipos<exp.size()-1; ipos++ ){
      exp[ipos-1] = stemp[ipos];
    }
    exp.resize( exp.size()-2 );
  }
  RemoveExpressionBracket(exp);
}



////////////////////////////////////////////////////////////////
// 引数expの中のもっとも優先順位の低いオペレターを検索する。
// itypeはオペレターの種類
// itype : 0 --> Numerical variable
// itype : 1 --> Argebric variable
// itype : 2 --> Unary Operator
// itype : 3 --> Binary Operator
int GetLowestPriorityOperator(
    unsigned int& ibegin,
    unsigned int& iend,
    int& itype,
    int& iopr,
    const std::string& exp )
{
  itype = -1;
  iopr = -1;
  
  const int ipriority_max = 10;
  
  unsigned int ipos_min = UINT_MAX;
  int ipriority_min = ipriority_max;
  int itype_min = -1;
  int iopr_min = -1;
  
  bool is_numeric = true;
  int iBracketDepth = 0;
  unsigned int iFirstBracket = UINT_MAX;
  for(unsigned int ipos_curr=0 ;ipos_curr<exp.size(); ipos_curr++ ){
    // 括弧内を読み飛ばす
    switch ( exp[ipos_curr] )
    {
      case '(':
        if( iFirstBracket == UINT_MAX ) iFirstBracket = ipos_curr;
        iBracketDepth++; continue;
      case ')':
        iBracketDepth--; continue;
      default:
        if( iBracketDepth != 0 ) continue;
    }
    assert( iBracketDepth == 0 );
    
    // ipriority_currを計算する
    int ipriority_curr;
    int itype_curr = -1;
    int iopr_curr = -1;
    switch( exp[ipos_curr] ){
      case '+':
        ibegin = ipos_curr;
        iend = ipos_curr+1;
        if( ipos_curr == 0 ){ // sign operator
          ipriority_curr=5;
          itype_curr=2;
          iopr_curr=0;
          break;
        }
        else if( exp[ipos_curr-1] != '*' && exp[ipos_curr-1] != '/' && exp[ipos_curr-1] != '^' ){
          iopr=0;
          itype=3;
          return 0;  // binary operator
        }
        continue;
      case '-':
        ibegin = ipos_curr;
        iend = ipos_curr+1;
        if( ipos_curr == 0 ){ // sign operator
          ipriority_curr=5;
          itype_curr=2;
          iopr_curr=1;
          break;
        }
        else if( exp[ipos_curr-1] != '*' && exp[ipos_curr-1] != '/' && exp[ipos_curr-1] != '^' ){
          iopr=1;
          itype=3;
          return 0;  // binary operator
        }
        continue;
      case '*':
        ipriority_curr=3;
        itype_curr=3;
        iopr_curr=2;
        break;
      case '/':
        ipriority_curr=4;
        itype_curr=3;
        iopr_curr=3;
        break;
      case '^':
        ipriority_curr=6;
        itype_curr=3;
        iopr_curr=4;
        break;
      default:
        ipriority_curr = ipriority_max;
        if( ( exp[ipos_curr] < '0' || exp[ipos_curr] > '9' ) && exp[ipos_curr] != '.' ){
          is_numeric = false;
        }
        break;
    }
    if ( ipriority_curr < ipriority_min ){
      ipos_min = ipos_curr;
      ipriority_min = ipriority_curr;
      itype_min = itype_curr;
      iopr_min = iopr_curr;
    }
  }
  
  if( iBracketDepth != 0 ){
    std::cout << "Error!-->Not corresponding bracket" << std::endl;
    return 2;
  }
  if( ( ipos_min==0 || ipos_min==exp.size()-1 ) && itype_min==3 ){
    std::cout << "Error!-->Binary operator misplaced " << std::endl;
    return 9;
  }
  if( iFirstBracket == 0 && ipos_min==-1 ){
    std::cout << "Error!-->(hoge)foo" << std::endl;
    return 3;
  }
  // binary operator(+とか-)
  if( ipos_min != UINT_MAX ){
//    std::cout << ipos_min << " " << itype_min << std::endl;
    assert( itype_min == 3 );
    assert( iopr_min != -1 );
    ibegin = ipos_min;
    iend = ipos_min+1;
    itype = itype_min;
    iopr = iopr_min;
    return 0;
  }
  // sign operator( -hoge )
  if( ipos_min == 0 ){
    assert( itype_min == 2 );
    assert( iopr_min != -1 );
    ibegin = ipos_min; iend = ipos_min+1;
    itype = itype_min; iopr = iopr_min;
    return 0;
  }
  // preposition operator ( sin, tanみたいの )
  if( !is_numeric && iFirstBracket != UINT_MAX ){
    if( exp[ exp.size()-1 ] != ')' ){  // avoid  sin(x+y)a
      std::cout << "Error!-->hoge1" << std::endl;
      return 1;
    }
    if( iFirstBracket == exp.size()-2 ){  // avoid sin()
      std::cout << "Error!-->hoge2" << std::endl;
      return 2;
    }
    std::string str1(exp,0,iFirstBracket);
    int iopr0 = CUnaryOperator::GetOprInd(str1);
    if( iopr0 != -1 ){
      ibegin = 0; iend = iFirstBracket;
      itype = 2; iopr = iopr0;
      return 0;
    }
    else{
      std::cout << "Error!-->Unregistered Unary Operator" << std::endl;
      return 1;
    }
  }
  // variable
  ibegin = 0;
  iend = static_cast<unsigned int>(exp.size());
  if( is_numeric ){ itype=0; iopr=-1; }  // 数値(5.4321みたいの)
  else{
    int iopr0 = COperand::GetOprInd(exp);
    if( iopr0 != -1 ){ itype = 1; iopr = iopr0; }  // PIみたいなの
    else{ itype=1; iopr=-1; }
  }
  return 0;
}

struct SExpCompo{
  std::string sOpe;
  int iOpeType;
  int iOpe;
};

bool MakeRPN(
    unsigned int icur_old,
    std::vector<SExpCompo>& exp_node_vec)
{
  assert( icur_old < exp_node_vec.size() );
  
  unsigned int ibegin0,iend0;
  int itype0,iopr0;
  {
    std::string& cur_exp = exp_node_vec[icur_old].sOpe;
    unsigned int itmp0, itmp1;
    int itmp2, itmp3;
    int ierr = GetLowestPriorityOperator(itmp0, itmp1, itmp2, itmp3, cur_exp);
    if( ierr != 0 ){
      std::cout << "Error!-->Cannot interprit this expression : " << cur_exp << std::endl;
      return false;
    }
    if( itmp0 == 0 && itmp1 == (int)cur_exp.size() ){  // case operand 算術数値(PIみたいなの)
      assert( itmp2 == 0 || itmp2 == 1 );
      assert( itmp3>=-1 && itmp3<=COperand::GetMaxOprInd() );
      exp_node_vec[icur_old].iOpeType = itmp2;
      exp_node_vec[icur_old].iOpe = itmp3;
      return true;
    }
    assert( itmp2 != 0 && itmp2 != 1 );
    ibegin0 = static_cast<unsigned int>(itmp0);
    iend0 = static_cast<unsigned int>(itmp1);
    assert( itmp2 >=0 && itmp2 <= 3 );
    itype0 = itmp2;
    assert( ( itmp2==2 && itmp3>=0 && itmp3<=CUnaryOperator::MaxOprInd() ) ||
           ( itmp2==3 && itmp3>=0 && itmp3<=CBinaryOperator::MaxOprInd() ) );
    iopr0 = itmp3;
  }
  
  std::string cur_old_exp = exp_node_vec[icur_old].sOpe;
  
  {
    unsigned int ileft = icur_old;
    SExpCompo& left_compo = exp_node_vec[ileft];
    left_compo.sOpe.assign( cur_old_exp, iend0, cur_old_exp.size()-iend0 );
    RemoveExpressionBracket(left_compo.sOpe);
    if( left_compo.sOpe.empty() ){ return false; }
    if( !MakeRPN(ileft,exp_node_vec) ){ return false; }
  }
  
  if( ibegin0 > 0 ){
    const unsigned int iright = static_cast<unsigned int>(exp_node_vec.size());
    exp_node_vec.resize( exp_node_vec.size()+1 );
    SExpCompo& right_compo = exp_node_vec[iright];
    right_compo.sOpe.assign( cur_old_exp, 0, ibegin0 );
    RemoveExpressionBracket(right_compo.sOpe);
    if( right_compo.sOpe.empty() ){ return false; }
    if( !MakeRPN(iright,exp_node_vec) ){ return false; }
  }
  
  {
    const size_t icur_new = exp_node_vec.size();
    exp_node_vec.resize( exp_node_vec.size()+1 );
    SExpCompo& new_compo = exp_node_vec[icur_new];
    new_compo.sOpe.assign( cur_old_exp, ibegin0, iend0-ibegin0 );
    new_compo.iOpeType = itype0;
    new_compo.iOpe = iopr0;
  }
  
  return true;
}


bool MakeCmdAry(
    std::vector<delfem2::evalmathexp::CCmd*>& cmd_vec,
    std::vector<dfm2::evalmathexp::CKey>& m_aKey,
    const std::vector<SExpCompo>& exp_node_vec)
{
  cmd_vec.resize( exp_node_vec.size(), nullptr );
  for(unsigned int iexp=0;iexp<exp_node_vec.size();iexp++){
    const SExpCompo& compo = exp_node_vec[iexp];
    if( compo.iOpeType == 0 ){ // numeric
      char* e;
      double val = strtod(compo.sOpe.c_str(),&e);
      if (errno != ERANGE){
        cmd_vec[iexp] = new COperand(val);
      }
      else if (val == HUGE_VAL){
        std::cout << "Exceeding the range of (double)" << std::endl;
        return false;
      }
    }
    else if( compo.iOpeType == 1 ){ // symbol
      if( compo.iOpe == -1 ){
        bool bflg = false;
        for(auto & ikey : m_aKey){
          if( ikey.m_Name == compo.sOpe ){
            ikey.m_aiCmd.push_back(iexp);
            bflg = true;
            break;
          }
        }
        if( !bflg ){
          std::cout << "次のオペランドを解釈できませんでした:" << compo.sOpe << std::endl;
          return false;
        }
        cmd_vec[iexp] = new COperand(0.0);
      }
      else{
        int iopr0 = compo.iOpe;
        if( iopr0 >=0 && iopr0 <= COperand::GetMaxOprInd() ){
          cmd_vec[iexp] = new COperand( COperand::GetValue(compo.iOpe) );
        }
        else{
          std::cout << "Not assumed Operand" << std::endl;
          assert(0);
        }
      }
    }
    else if( compo.iOpeType == 2 ){ // unary
      cmd_vec[iexp] = new CUnaryOperator( compo.iOpe );
    }
    else if( compo.iOpeType == 3 ){ // binary
      cmd_vec[iexp] = new CBinaryOperator( compo.iOpe );
    }
    else{
      std::cout << "Error!--> " << compo.sOpe << " " << compo.iOpeType << std::endl;
      assert(0);
    }
  }
  return true;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

dfm2::CMathExpressionEvaluator::~CMathExpressionEvaluator(){
	for(auto & icmd : m_apCmd){
		delete icmd;
	}
}

void dfm2::CMathExpressionEvaluator::SetKey
 (const std::string& key_name, double key_val)
{
	for(auto & ikey : m_aKey){
		if( ikey.m_Name == key_name){
			ikey.m_Val = key_val;
			for(unsigned int icmd0 : ikey.m_aiCmd){	// このKeyをもっている全てのCmdに値をセット
			  assert( icmd0 < m_apCmd.size() );
				m_apCmd[icmd0]->SetValue( key_val );
			}
			return;
		}
	}
	// 未登録のKeyの場合は以下のルーティン
	m_aKey.emplace_back(key_name,key_val );
}

bool dfm2::CMathExpressionEvaluator::SetExp(const std::string& exp){
  m_is_valid = false;
	m_sExp = exp;
	// 消去
	for(auto & icmd : m_apCmd){ delete icmd; }
	for(auto & ikey : m_aKey){ ikey.m_aiCmd.clear(); }
	m_apCmd.clear();

	////////////////
	std::string tmp_exp = exp;
	RemoveExpressionSpaceNewline(tmp_exp);
	RemoveExpressionBracket(tmp_exp);
	if( tmp_exp.empty() ){
		m_is_valid = false;
		return false;
	}
	{
		std::vector<SExpCompo> exp_vec;
		exp_vec.reserve(exp.size());
		exp_vec.resize(1);
		exp_vec[0].sOpe = tmp_exp;
		if( !MakeRPN(0,exp_vec) ){
			m_is_valid = false;
			return false;
		}

/*		{	// 逆ポーランド表記を表示
			for(unsigned int iexno=0;iexno<exp_vec.size();iexno++){
				std::cout << exp_vec[iexno].sOpe << " ";
			}
			std::cout << std::endl;
		}*/

    if( !MakeCmdAry(m_apCmd, m_aKey,
                    exp_vec) ){
			m_is_valid = false;
			return false;
		}
	
/*
		int iflg = true;
		int stack_size = 0;
		for(int icmd=0;icmd<m_apCmd.size();icmd++){
			if( dynamic_cast<COperand*>(m_apCmd[icmd]) != NULL ){
				stack_size++;
			}
			else if( dynamic_cast<CUnaryOperator*>(m_apCmd[icmd]) != NULL ){
				if( stack_size <= 0 ) iflg = false; break;
			}
			else if( dynamic_cast<CBinaryOperator*>(m_apCmd[icmd]) != NULL ){
				if( stack_size <= 1 ) iflg = false;	break;
			}
		}
		if( !iflg ){
			std::cout << "Error!-->Operator and Operand mismatch" << std::endl;
			m_bValid = false;
			return false;
		}
*/	
	}
	{	// 登録されているキー全ての値を設定
		for(const auto & key : m_aKey){
		  double val = key.m_Val;
			for(unsigned int icmd0 : key.m_aiCmd){	// このKeyをもっている全てのCmdに値をセット
					assert( icmd0 < m_apCmd.size() );
				m_apCmd[icmd0]->SetValue( val );
			}
		}	
	}
	m_is_valid = true;
	return true;
}

double dfm2::CMathExpressionEvaluator::Eval() const{ // evaluating the math expression
	std::vector<double> stack;
	stack.reserve(128);
	for(auto icmd : m_apCmd){
		icmd->DoOperation(stack);
	}
	assert( stack.size() == 1 );
	return stack[0];
}
