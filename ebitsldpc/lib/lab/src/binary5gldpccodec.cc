#include "binary5gldpccodec.h"

namespace lab {
Binary5GLDPCCodec::Binary5GLDPCCodec(const Binary5GLDPCCodec &codec)
    : BinaryLDPCCodec(codec), code_len_no_puncture_(codec.code_len_no_puncture_),
      code_len_puncture_(codec.code_len_puncture_), lifting_factor_(codec.lifting_factor_) {
  cc_no_puncture_ = new int[code_len_no_puncture_];
  cc_soft_no_puncture_ = new double[code_len_no_puncture_];
}

Binary5GLDPCCodec::Binary5GLDPCCodec(const toml::value &arguments)
    : code_len_no_puncture_(0), code_len_puncture_(0), lifting_factor_(0),
      cc_no_puncture_(nullptr), cc_soft_no_puncture_(nullptr) {
  int i, j;
  int row_no, row_deg, col_no;
  char temp_str[80] = {' '};
  FILE *fp;
  const auto ldpc = toml::find(arguments, "ldpc");
  max_iter_ = toml::find<int>(ldpc, "max_iter");
  encoder_active_ = toml::find<bool>(ldpc, "active");
  std::string matrix_file = toml::find<std::string>(ldpc, "matrix_file");

  //Read H from file temp_str
  if ((fp = fopen(matrix_file.c_str(), "r")) == nullptr) {
    lab::logger::ERROR("Cannot Opne" + matrix_file, true);
    exit(-1);
  }
  fscanf(fp, "%s", temp_str);
  fscanf(fp, "%d %d %d %d", &num_row_, &num_col_, &code_chk_, &lifting_factor_);
  code_len_no_puncture_ = num_col_;
  code_len_puncture_ = num_col_ - lifting_factor_ * 2;
  code_dim_ = code_len_no_puncture_ - code_chk_;
  coderate_ = (double)code_dim_ / code_len_puncture_;
  cc_no_puncture_ = new int[code_len_no_puncture_];
  cc_soft_no_puncture_ = new double[code_len_no_puncture_];
  row_head_ = new Edge[num_row_];
  col_head_ = new Edge[num_col_];
  syndrom_soft_ = new double[num_row_];
  for (i = 0; i < num_row_; i++) {
    (row_head_ + i)->m_row_no = i;
    (row_head_ + i)->m_col_no = -1;
    (row_head_ + i)->left = row_head_ + i;
    (row_head_ + i)->right = row_head_ + i;
    (row_head_ + i)->up = row_head_ + i;
    (row_head_ + i)->down = row_head_ + i;
  }
  for (i = 0; i < num_col_; i++) {
    (col_head_ + i)->m_row_no = -1;
    (col_head_ + i)->m_col_no = i;
    (col_head_ + i)->left = col_head_ + i;
    (col_head_ + i)->right = col_head_ + i;
    (col_head_ + i)->up = col_head_ + i;
    (col_head_ + i)->down = col_head_ + i;
  }
  fscanf(fp, "%s", temp_str);
  Edge *temp_edge;
  for (i = 0; i < num_row_; i++) {
    fscanf(fp, "%d %d", &row_no, &row_deg);
    for (j = 0; j < row_deg; j++) {
      temp_edge = new Edge;
      temp_edge->m_row_no = row_no;
      fscanf(fp, "%d", &col_no);
      temp_edge->m_col_no = col_no;
      temp_edge->right = (row_head_ + i)->right;
      (row_head_ + i)->right = temp_edge;
      temp_edge->left = row_head_ + i;
      (temp_edge->right)->left = temp_edge;
      temp_edge->down = (col_head_ + col_no)->down;
      (col_head_ + col_no)->down = temp_edge;
      temp_edge->up = col_head_ + col_no;
      (temp_edge->down)->up = temp_edge;
    }
  }
  fclose(fp);
  if (encoder_active_ == 1)
    SystemMatrixH();
  cc_hat_ = new int[code_len_no_puncture_];
}

Binary5GLDPCCodec::~Binary5GLDPCCodec() {
  delete[] cc_no_puncture_;
  delete[] cc_soft_no_puncture_;
}

void
Binary5GLDPCCodec::Encoder(int *uu, int *cc) const {
  int i, j, t;
  for (i = 0; i < code_len_puncture_; i++)
    cc_no_puncture_[i] = 0;

  //codeword = [parity_check_bits information_bits]
  if (encoder_active_) {
    for (t = 0; t < code_dim_; t++)
      cc_no_puncture_[t] = uu[t];
    for (t = 0; t < code_chk_; t++) {
      cc_no_puncture_[code_dim_ + t] = 0;
      for (j = 0; j < code_dim_; j++)
        cc_no_puncture_[code_dim_ + t] ^= (cc_no_puncture_[j] & enc_h_[t][j]);
    }
    for (t = 0; t < code_len_puncture_; t++) {
      cc[t] = cc_no_puncture_[t + lifting_factor_ * 2];
    }
  } else {
    for (i = 0; i < code_dim_; i++)
      uu[i] = 0;
    for (i = 0; i < code_len_puncture_; i++)
      cc[i] = 0;
  }
}

int
Binary5GLDPCCodec::Decoder(const double *M2V, int *uu_hat, int iter_cnt) {
  int i;
  int iter;
  int parity_check;
  double temp0, temp1, temp_sum;
  Edge *p_edge;

  //initialization
  InitMsg();
  //iteration
  for (iter = 0; iter < iter_cnt; iter++) {
    //from vnode to cnode
    for (i = 0; i < num_col_; i++) {
      //forward
      if (i < lifting_factor_ * 2) {
        p_edge = (col_head_ + i)->down;
        p_edge->m_alpha[0] = 0.5;
        p_edge->m_alpha[1] = 1.0 - 0.5;
      } else {
        p_edge = (col_head_ + i)->down;
        p_edge->m_alpha[0] = M2V[i - lifting_factor_ * 2];
        p_edge->m_alpha[1] = 1.0 - M2V[i - lifting_factor_ * 2];
      }
      while (p_edge->m_row_no != -1) {
        p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
        p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
        temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
        p_edge->down->m_alpha[0] /= temp_sum;
        p_edge->down->m_alpha[1] /= temp_sum;
        p_edge = p_edge->down;
      }
      //hard decision
      if ((col_head_ + i)->m_alpha[0] > (col_head_ + i)->m_alpha[1])
        cc_hat_[i] = 0;
      else
        cc_hat_[i] = 1;

      //backward
      p_edge = (col_head_ + i)->up;
      p_edge->m_beta[0] = 1.0;
      p_edge->m_beta[1] = 1.0;
      while (p_edge->m_row_no != -1) {
        temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0];
        temp1 = p_edge->m_alpha[1] * p_edge->m_beta[1];
        temp_sum = temp0 + temp1;
        p_edge->m_v2c[0] = temp0 / temp_sum;
        p_edge->m_v2c[1] = temp1 / temp_sum;
        p_edge->up->m_beta[0] = p_edge->m_beta[0] * p_edge->m_c2v[0];
        p_edge->up->m_beta[1] = p_edge->m_beta[1] * p_edge->m_c2v[1];
        temp_sum = p_edge->up->m_beta[0] + p_edge->up->m_beta[1];
        p_edge->up->m_beta[0] /= temp_sum;
        p_edge->up->m_beta[1] /= temp_sum;
        p_edge = p_edge->up;
      }
    }
    for (i = 0; i < code_dim_; i++) {
      //			uu_hat[i] = cc_hat_[i+m_codechk];
      uu_hat[i] = cc_hat_[i];//kite使用
    }
    //parity checking
    success_ = 1;
    for (i = 0; i < num_row_; i++) {
      parity_check = 0;
      p_edge = (row_head_ + i)->right;
      while (p_edge->m_col_no != -1) {
        parity_check = parity_check ^ cc_hat_[p_edge->m_col_no];
        p_edge = p_edge->right;
      }
      if (parity_check != 0) {
        success_ = 0;
        break;
      }
    }
    if (success_ == 1)
      break;

    //from c node to v node
    for (i = 0; i < num_row_; i++) {
      //forward
      p_edge = (row_head_ + i)->right;
      p_edge->m_alpha[0] = 1.0;
      p_edge->m_alpha[1] = 0.0;
      while (p_edge->m_col_no != -1) {//over trellis with two states
        p_edge->right->m_alpha[0] =
            p_edge->m_alpha[0] * p_edge->m_v2c[0] + p_edge->m_alpha[1] * p_edge->m_v2c[1];
        p_edge->right->m_alpha[1] =
            p_edge->m_alpha[0] * p_edge->m_v2c[1] + p_edge->m_alpha[1] * p_edge->m_v2c[0];
        temp_sum = p_edge->right->m_alpha[0] + p_edge->right->m_alpha[1];
        p_edge->right->m_alpha[0] /= temp_sum;
        p_edge->right->m_alpha[1] /= temp_sum;
        p_edge = p_edge->right;
      }
      //backward
      p_edge = (row_head_ + i)->left;
      p_edge->m_beta[0] = 1.0;
      p_edge->m_beta[1] = 0.0;
      while (p_edge->m_col_no != -1) {
        temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0] + p_edge->m_alpha[1] * p_edge->m_beta[1];
        temp1 = p_edge->m_alpha[0] * p_edge->m_beta[1] + p_edge->m_alpha[1] * p_edge->m_beta[0];
        temp_sum = temp0 + temp1;
        p_edge->m_c2v[0] = temp0 / temp_sum;
        p_edge->m_c2v[1] = temp1 / temp_sum;
        if (p_edge->m_c2v[0] > 1.0 - kSmallestProb)
          p_edge->m_c2v[0] = 1.0 - kSmallestProb;
        if (p_edge->m_c2v[0] < kSmallestProb)
          p_edge->m_c2v[0] = kSmallestProb;
        p_edge->m_c2v[1] = 1.0 - p_edge->m_c2v[0];
        p_edge->left->m_beta[0] =
            p_edge->m_beta[0] * p_edge->m_v2c[0] + p_edge->m_beta[1] * p_edge->m_v2c[1];
        p_edge->left->m_beta[1] =
            p_edge->m_beta[0] * p_edge->m_v2c[1] + p_edge->m_beta[1] * p_edge->m_v2c[0];
        temp_sum = p_edge->left->m_beta[0] + p_edge->left->m_beta[1];
        p_edge->left->m_beta[0] /= temp_sum;
        p_edge->left->m_beta[1] /= temp_sum;
        p_edge = p_edge->left;
      }
      syndrom_soft_[i] = (row_head_ + i)->m_alpha[0];
    }
  }
  return iter + (iter < max_iter_);
}

int
Binary5GLDPCCodec::code_len_puncture() const {
  return code_len_puncture_;
}

void
Binary5GLDPCCodec::SystemMatrixH() {
  int i, j, ii, jj, m, n;
  int temp;
  int flag;
  int *tempP;
  char **tempH;
  Edge *p_edge;
  dec_h_ = new char *[num_row_];
  for (i = 0; i < num_row_; i++)
    dec_h_[i] = new char[num_col_];
  for (i = 0; i < num_row_; i++) {
    for (j = 0; j < num_col_; j++)
      dec_h_[i][j] = 0;
  }
  for (i = 0; i < num_row_; i++) {
    p_edge = (row_head_ + i)->right;
    while (p_edge->m_col_no != -1) {
      dec_h_[i][p_edge->m_col_no] = 1;
      p_edge = p_edge->right;
    }
  }
  tempP = new int[num_col_];
  for (j = 0; j < num_col_; j++)
    tempP[j] = j;
  tempH = new char *[num_row_];
  for (i = 0; i < num_row_; i++)
    tempH[i] = new char[num_col_];
  for (i = 0; i < num_row_; i++) {
    for (j = 0; j < num_col_; j++)
      tempH[i][j] = dec_h_[i][j];
  }
  enc_h_ = new char *[num_row_];
  for (i = 0; i < num_row_; i++) {
    enc_h_[i] = new char[num_col_];
  }
  for (i = 0; i < num_row_; i++) {
    for (j = 0; j < num_col_; j++) {
      enc_h_[i][j] = dec_h_[i][j];
    }
  }
  code_chk_ = 0;
  for (i = num_row_ - 1; i >= 0; --i) {
    flag = 0;
    for (jj = i + num_col_ - num_row_; jj >= 0; --jj) {
      for (ii = i; ii >= 0; --ii) {
        if (enc_h_[ii][jj] != 0) {
          flag = 1;
          break;
        }
      }
      if (flag == 1) {
        code_chk_++;
        break;
      }
    }
    if (flag == 0)
      break;
    else {
      //swap i and ii row
      if (ii != i) {
        for (n = 0; n < num_col_; n++) {
          temp = enc_h_[i][n];
          enc_h_[i][n] = enc_h_[ii][n];
          enc_h_[ii][n] = temp;
        }
      }
      //swap (i + m_num_col - m_num_row) and jj col
      if (jj != i + num_col_ - num_row_) {
        temp = tempP[i + num_col_ - num_row_];
        tempP[i + num_col_ - num_row_] = tempP[jj];
        tempP[jj] = temp;
        for (m = 0; m < num_row_; m++) {
          temp = enc_h_[m][i + num_col_ - num_row_];
          enc_h_[m][i + num_col_ - num_row_] = enc_h_[m][jj];
          enc_h_[m][jj] = temp;
        }
      }
      //elimination
      for (m = num_row_ - 1; m >= 0; --m) {
        if (m != i && enc_h_[m][i + num_col_ - num_row_] == 1) {
          for (n = 0; n < num_col_; n++)
            enc_h_[m][n] ^= enc_h_[i][n];
        }
      }
    }
  }

  //for (i = 0; i < m_num_col; i++)
  //{
  //	if (tempP[i] != i)
  //	{
  //		fprintf(stderr, "\nWarning: please make sure the information bits on the right!");
  //	}
  //}

  for (j = 0; j < num_col_; j++) {
    for (i = 0; i < num_row_; i++) {
      dec_h_[i][j] = tempH[i][tempP[j]];
    }
  }

  /* By Zhaosc
  -----由于在Gauss消元过程中有行列变换
  -----需要删除原有的Tanner Graph，重新建立.
  */
  FreeTannerGraph();
  row_head_ = new Edge[num_row_];
  col_head_ = new Edge[num_col_];
  for (i = 0; i < num_row_; i++) {
    (row_head_ + i)->m_row_no = i;
    (row_head_ + i)->m_col_no = -1;
    (row_head_ + i)->left = row_head_ + i;
    (row_head_ + i)->right = row_head_ + i;
    (row_head_ + i)->up = row_head_ + i;
    (row_head_ + i)->down = row_head_ + i;
  }
  for (i = 0; i < num_col_; i++) {
    (col_head_ + i)->m_row_no = -1;
    (col_head_ + i)->m_col_no = i;
    (col_head_ + i)->left = col_head_ + i;
    (col_head_ + i)->right = col_head_ + i;
    (col_head_ + i)->up = col_head_ + i;
    (col_head_ + i)->down = col_head_ + i;
  }
  for (i = 0; i < num_row_; i++) {
    for (j = 0; j < num_col_; j++) {
      if (dec_h_[i][j] != 0) {
        p_edge = new Edge;
        p_edge->m_row_no = i;
        p_edge->m_col_no = j;
        p_edge->right = (row_head_ + i)->right;
        (row_head_ + i)->right = p_edge;
        p_edge->left = row_head_ + i;
        (p_edge->right)->left = p_edge;
        p_edge->down = (col_head_ + j)->down;
        (col_head_ + j)->down = p_edge;
        p_edge->up = col_head_ + j;
        (p_edge->down)->up = p_edge;
      }
    }
  }
  code_len_no_puncture_ = num_col_;
  code_dim_ = code_len_no_puncture_ - code_chk_;
  coderate_ = (double)code_dim_ / code_len_puncture_;
  delete[] tempP;
  for (i = 0; i < num_row_; i++) {
    delete[] tempH[i];
    delete[] dec_h_[i];
  }
  delete[] tempH;
  delete[] dec_h_;
}
}// namespace lab
